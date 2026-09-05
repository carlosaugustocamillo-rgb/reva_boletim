import json
import tempfile
import unittest
from datetime import date
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import requests

from pubmed_related import (
    PubMedRelatedClient, RelatedConfig, article_keys, context_for_script, date_bounds,
    enabled_from_env, enrich_episode, parse_articles, prefilter, publication_date,
    references_for_notes, save_json, screen_candidates,
)


ABSTRACT = "Adults with breast cancer received supervised exercise during chemotherapy. Fatigue did not differ between groups."
ANCHOR = {
    "pmid": "100", "titulo": "A new exercise trial", "doi": "10.123/new",
    "data_publicacao": "2026-09-04", "resumo_original": ABSTRACT,
    "resumo_traduzido": "Resumo traduzido", "autores": ["Silva A"], "tipos": ["Randomized Controlled Trial"],
}


def candidate(pmid="200", **changes):
    return {"pmid": pmid, "titulo": f"Earlier exercise trial {pmid}", "doi": f"10.123/{pmid}",
            "data_publicacao": "2025-01-01", "resumo_original": ABSTRACT, "autores": ["Jones B"],
            "tipos": ["Randomized Controlled Trial"], "similarity_score": 50, "mesh": ["Humans"], **changes}


def decision(pmid="200", **changes):
    return {"pmid": pmid, "incluir": True, "populacao": "compativel", "intervencao": "compativel",
            "desfecho": "compativel", "comparador": "compativel", "desenho": "ensaio_randomizado",
            "justificativa": "Mesma população e intervenção; resultado nulo preservado.",
            "trecho_ancora": "Adults with breast cancer received supervised exercise during chemotherapy.",
            "trecho_candidato": "Adults with breast cancer received supervised exercise during chemotherapy.",
            "diferenca_importante": "Diferentes doses de exercício.", **changes}


def llm(decisions, finish_reason="stop"):
    client = MagicMock()
    client.with_options.return_value.chat.completions.create.return_value = SimpleNamespace(
        choices=[SimpleNamespace(finish_reason=finish_reason, message=SimpleNamespace(content=json.dumps({"decisoes": decisions})))],
        usage=SimpleNamespace(prompt_tokens=150, completion_tokens=80),
    )
    return client


XML = b'''<PubmedArticleSet><PubmedArticle><MedlineCitation><PMID>200</PMID><Article>
<ArticleTitle>Exercise <i>during</i> treatment</ArticleTitle><Journal><Title>Journal</Title><JournalIssue>
<PubDate><Year>2025</Year><Month>Feb</Month></PubDate></JournalIssue></Journal>
<ArticleDate DateType="Electronic"><Year>2025</Year><Month>01</Month><Day>15</Day></ArticleDate>
<Abstract><AbstractText Label="METHODS">Adults with breast cancer received <b>exercise</b>.</AbstractText>
<AbstractText Label="RESULTS">No difference in fatigue.</AbstractText></Abstract>
<AuthorList><Author><CollectiveName>Exercise Group</CollectiveName></Author></AuthorList>
<PublicationTypeList><PublicationType>Randomized Controlled Trial</PublicationType></PublicationTypeList>
</Article><MeshHeadingList><MeshHeading><DescriptorName>Humans</DescriptorName></MeshHeading></MeshHeadingList>
</MedlineCitation><PubmedData><ArticleIdList><ArticleId IdType="doi">https://doi.org/10.123/ABC</ArticleId>
</ArticleIdList></PubmedData></PubmedArticle></PubmedArticleSet>'''


def links(ids=("200",), error=""):
    payload = "".join(f"<Link><Id>{pmid}</Id><Score>{index}</Score></Link>" for index, pmid in enumerate(ids))
    return f'<eLinkResult>{error}<LinkSet><IdList><Id>100</Id></IdList><LinkSetDb><LinkName>pubmed_pubmed</LinkName>{payload}</LinkSetDb></LinkSet></eLinkResult>'.encode()


def response(content=XML, status=200, headers=None):
    return SimpleNamespace(content=content, status_code=status, headers=headers or {})


class PubmedRelatedTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.sleep_patch = patch("pubmed_related.time.sleep")
        self.sleep = self.sleep_patch.start()
        self.addCleanup(self.sleep_patch.stop)
        self.session = MagicMock()
        self.pubmed = PubMedRelatedClient(self.temp.name, session=self.session, email="test@example.org")

    def test_parses_structured_abstract_inline_markup_author_and_doi(self):
        article = parse_articles(XML)[0]
        self.assertEqual(article["titulo"], "Exercise during treatment")
        self.assertEqual(article["doi"], "10.123/abc")
        self.assertEqual(article["autores"], ["Exercise Group"])
        self.assertEqual(article["data_publicacao"], "2025-01-15")
        self.assertEqual(article["resumo_original"], "METHODS: Adults with breast cancer received exercise.\n\nRESULTS: No difference in fatigue.")

    def test_partial_dates_and_ranges_are_not_invented(self):
        self.assertEqual(date_bounds("2024-02")[1], date(2024, 2, 29))
        self.assertEqual(date_bounds("2025 Dec-2026 Jan"), (date(2025, 1, 1), date(2026, 12, 31)))
        self.assertIsNone(date_bounds("unknown"))
        self.assertIsNone(date_bounds("2026-99-99"))
        self.assertEqual(publication_date({"Journal": {"JournalIssue": {"PubDate": {"Year": "2025"}}}}), "2025")

    def test_neighbors_batch_and_cache(self):
        self.session.get.side_effect = [response(links(("100", "200", "200"))), response()]
        items, cached = self.pubmed.candidates("100")
        self.assertFalse(cached)
        self.assertEqual([item["pmid"] for item in items], ["200"])
        calls = self.session.get.call_args_list
        self.assertEqual(calls[0].kwargs["params"]["cmd"], "neighbor_score")
        self.assertEqual(calls[0].kwargs["params"]["linkname"], "pubmed_pubmed")
        self.assertEqual(calls[1].kwargs["params"]["id"], "200")
        self.assertNotIn("api_key", calls[0].kwargs["params"])
        again, cached = self.pubmed.candidates("100")
        self.assertTrue(cached)
        self.assertEqual(again, items)
        self.assertEqual(self.session.get.call_count, 2)

    def test_sorts_neighbors_by_score_not_response_order(self):
        self.session.get.side_effect = [response(links(("200", "201"))), response(XML.replace(b"200", b"201").replace(b"</PubmedArticleSet>", XML.split(b"<PubmedArticleSet>")[1]))]
        items, _ = self.pubmed.candidates("100")
        self.assertEqual([item["pmid"] for item in items], ["201", "200"])

    def test_ignores_other_link_types(self):
        self.session.get.return_value = response(links().replace(b"pubmed_pubmed", b"pubmed_pubmed_refs"))
        self.assertEqual(self.pubmed.candidates("100")[0], [])
        self.assertEqual(self.session.get.call_count, 1)

    def test_empty_results_are_cached_briefly(self):
        self.session.get.return_value = response(links(()))
        self.assertEqual(self.pubmed.candidates("100")[0], [])
        self.assertTrue(self.pubmed.candidates("100")[1])
        path = Path(self.temp.name) / "v1_100_12.json"
        data = json.loads(path.read_text())
        data["saved_at"] -= 3601
        save_json(path, data)
        self.assertFalse(self.pubmed.candidates("100")[1])

    def test_corrupt_and_expired_cache_are_refetched(self):
        path = Path(self.temp.name) / "v1_100_12.json"
        save_json(path, {"saved_at": 0, "articles": [candidate()]})
        self.session.get.return_value = response(links(()))
        self.assertFalse(self.pubmed.candidates("100")[1])
        save_json(path, ["bad schema"])
        self.assertFalse(self.pubmed.candidates("100")[1])

    def test_rate_limit_backoff_and_timeout(self):
        self.session.get.side_effect = [response(status=429, headers={"Retry-After": "3"}), response(links(()))]
        self.pubmed.candidates("100")
        self.sleep.assert_any_call(3)
        self.assertEqual(self.session.get.call_count, 2)
        self.assertEqual(self.session.get.call_args.kwargs["timeout"], (5, 10))

    def test_timeout_retries_are_bounded_and_do_not_leak_secrets(self):
        self.session.get.side_effect = requests.Timeout("https://example?api_key=SECRET")
        with self.assertRaisesRegex(RuntimeError, "três tentativas") as error:
            self.pubmed.candidates("100")
        self.assertNotIn("SECRET", str(error.exception))
        self.assertEqual(self.session.get.call_count, 3)

    def test_long_retry_after_skips_optional_step(self):
        self.session.get.return_value = response(status=429, headers={"Retry-After": "120"})
        with self.assertRaises(RuntimeError):
            self.pubmed.candidates("100")
        self.assertEqual(self.session.get.call_count, 1)

    def test_xml_error_not_cached(self):
        self.session.get.return_value = response(links(error="<ERROR>temporarily unavailable</ERROR>"))
        with self.assertRaises(ValueError):
            self.pubmed.candidates("100")
        self.assertFalse(list(Path(self.temp.name).glob("*.json")))

    def test_bad_pmid_never_becomes_cache_path_or_request(self):
        with self.assertRaises(ValueError):
            self.pubmed.candidates("../secret")
        self.session.get.assert_not_called()

    def test_filter_deduplicates_pmid_doi_and_title(self):
        candidates = [candidate(pmid="100"), candidate(doi="https://doi.org/10.123/NEW"),
                      candidate(titulo=ANCHOR["titulo"].upper()), candidate(), candidate("201", doi="10.123/200")]
        valid, rejected = prefilter(candidates, ANCHOR, set(), date(2026, 9, 4))
        self.assertEqual([item["pmid"] for item in valid], ["200"])
        self.assertEqual(len(rejected), 4)

    def test_title_deduplication_ignores_accents_and_punctuation(self):
        anchor = {**ANCHOR, "titulo": "Exercício e câncer: reabilitação"}
        duplicate = candidate(titulo="Exercicio e cancer - reabilitacao")
        self.assertTrue(article_keys(anchor) & article_keys(duplicate))

    def test_filters_missing_abstract_dates_retractions_animals(self):
        changes = [dict(resumo_original=""), dict(data_publicacao=""), dict(data_publicacao="2027"),
                   dict(data_publicacao="2026"), dict(data_publicacao="2026-09-04"),
                   dict(tipos=["Retracted Publication"]), dict(alertas=["RetractionIn"]),
                   dict(mesh=["Animals"]), dict(tipos=["Preprint"])]
        for fields in changes:
            with self.subTest(fields=fields):
                valid, rejected = prefilter([candidate(**fields)], ANCHOR, set(), date(2026, 9, 4))
                self.assertEqual(valid, [])
                self.assertEqual(len(rejected), 1)

    def test_missing_anchor_date_never_fabricates_chronology(self):
        valid, rejected = prefilter([candidate()], {**ANCHOR, "data_publicacao": ""}, set(), date(2026, 9, 4))
        self.assertEqual(valid, [])
        self.assertEqual(rejected[0]["motivo"], "data_desconhecida")

    def test_screen_accepts_supported_negative_findings_and_preserves_metadata(self):
        client = llm([decision()])
        selected, audit, tokens = screen_candidates(ANCHOR, [candidate()], client, RelatedConfig())
        self.assertEqual(selected[0]["resumo_original"], ABSTRACT)
        self.assertTrue(audit[0]["aceito_validacao"])
        self.assertEqual(tokens["prompt_tokens"], 150)
        options = client.with_options.return_value.chat.completions.create.call_args.kwargs
        self.assertTrue(options["response_format"]["json_schema"]["strict"])
        self.assertEqual(options["temperature"], 0)

    def test_screen_rejects_uncertain_pico_unknown_design_or_fabricated_quote(self):
        for changes in [dict(populacao="incerta"), dict(intervencao="incompativel"),
                        dict(desenho="outro"), dict(trecho_candidato="A fabricated quotation about efficacy."),
                        dict(trecho_ancora="Adults"), dict(incluir=False)]:
            with self.subTest(changes=changes):
                selected, _, _ = screen_candidates(ANCHOR, [candidate()], llm([decision(**changes)]), RelatedConfig())
                self.assertEqual(selected, [])

    def test_screen_fails_closed_on_missing_duplicate_or_invented_ids(self):
        for decisions in [[], [decision("999")], [decision(), decision()]]:
            with self.subTest(decisions=decisions), self.assertRaises(ValueError):
                screen_candidates(ANCHOR, [candidate()], llm(decisions), RelatedConfig())

    def test_screen_fails_closed_on_truncated_model_output(self):
        with self.assertRaises(ValueError):
            screen_candidates(ANCHOR, [candidate()], llm([decision()], finish_reason="length"), RelatedConfig())

    def test_disabled_never_calls_services(self):
        client, pubmed = MagicMock(), MagicMock()
        result = enrich_episode([ANCHOR], self.temp.name, client, enabled=False, pubmed=pubmed)
        self.assertEqual(result["status"], "desativado")
        pubmed.candidates.assert_not_called()
        client.with_options.assert_not_called()
        self.assertEqual(references_for_notes(result), "")

    def test_episode_caps_anchors_and_deduplicates_across_them(self):
        pubmed = MagicMock()
        pubmed.candidates.return_value = ([candidate()], False)
        articles = [ANCHOR, {**ANCHOR, "pmid": "101", "titulo": "Second trial", "doi": "10.123/second"}, {**ANCHOR, "pmid": "102"}]
        report = enrich_episode(articles, self.temp.name, llm([decision()]), pubmed=pubmed, today=date(2026, 9, 4), log=lambda _: None)
        self.assertEqual(pubmed.candidates.call_count, 2)
        self.assertEqual(len(report["estudos"][0]["referencias"]), 1)
        self.assertEqual(report["estudos"][1]["referencias"], [])
        self.assertNotIn("contexto_pubmed", ANCHOR)
        notes = references_for_notes(report)
        self.assertIn("https://pubmed.ncbi.nlm.nih.gov/200/", notes)
        self.assertIn("Contexto anterior para PMID 100", notes)

    def test_failure_in_discovery_or_screening_keeps_episode_without_context(self):
        for discovery_fails in (True, False):
            pubmed = MagicMock()
            if discovery_fails:
                pubmed.candidates.side_effect = requests.Timeout("secret")
            else:
                pubmed.candidates.return_value = ([candidate()], False)
            client = MagicMock()
            client.with_options.side_effect = RuntimeError("secret")
            report = enrich_episode([ANCHOR], self.temp.name, client, pubmed=pubmed, log=lambda _: None)
            self.assertEqual(report["status"], "parcial")
            self.assertEqual(report["estudos"][0]["referencias"], [])
            self.assertNotIn("secret", json.dumps(report))

    def test_reference_limit_after_screening(self):
        articles = [candidate(), candidate("201", similarity_score=100), candidate("202", similarity_score=10)]
        selected, _, _ = screen_candidates(ANCHOR, articles, llm([decision(item["pmid"]) for item in articles]), RelatedConfig())
        self.assertEqual([item["pmid"] for item in selected], ["201", "200"])

    def test_context_is_grounded_and_no_context_does_not_claim_no_evidence(self):
        self.assertIn("NÃO demonstra inexistência", context_for_script(None))
        context = context_for_script({"referencias": [{**candidate(), "triagem": decision()}]})
        self.assertIn(ABSTRACT, context)
        self.assertIn("não novidades publicadas nesta semana", context)

    def test_configuration_bounds_and_global_off_switch(self):
        with patch.dict("os.environ", {"PODCAST_PUBMED_RELATED_MAX_ANCHORS": "100", "PODCAST_PUBMED_RELATED_MAX_CANDIDATES": "bad", "PODCAST_PUBMED_RELATED_ENABLED": "false"}):
            self.assertEqual(RelatedConfig.from_env().max_anchors, 6)
            self.assertEqual(RelatedConfig.from_env().max_candidates, 12)
            self.assertFalse(enabled_from_env())


if __name__ == "__main__":
    unittest.main()
