"""Opt-in paid text comparison. Never imports the pipeline or calls TTS/publishing.

Run with the project's Python: scripts/compare_podcast_conversation.py
    --before-prompts /path/to/snapshot --output /private/tmp/comparison
Both versions use identical sources, word target, model and generation settings.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import statistics
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import podcast_editorial as editorial


def rhythm(script):
    """Descriptive counters only: these are not a score or a naturalness gate."""
    turns = [turn for block in editorial.dialogues({"script": script}) for turn in block]
    sentences = [sentence for turn in turns
                 for sentence in re.split(r"(?<=[.!?])\s+", turn["text"].strip()) if sentence]
    lengths = [len(sentence.split()) for sentence in sentences]
    return {"words": sum(len(turn["text"].split()) for turn in turns),
            "turns": len(turns), "mean_words_per_sentence": round(statistics.mean(lengths), 1),
            "sentences_over_25_words": sum(size > 25 for size in lengths),
            "max_sentence_words": max(lengths),
            "note": "Contagens aproximadas; requerem avaliação humana da interação e dos fatos."}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--before-prompts", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--only", choices=("before", "after"),
                        help="Roda só uma versão, por exemplo ao refinar o prompt após a primeira comparação.")
    args = parser.parse_args()
    if args.output.exists():
        parser.error("Use uma pasta nova para preservar comparações anteriores.")
    for stage in ("plan", "write", "audit"):
        if not (args.before_prompts / f"{stage}.txt").is_file():
            parser.error(f"Snapshot anterior sem {stage}.txt")

    from dotenv import load_dotenv
    from openai import OpenAI
    load_dotenv(ROOT / ".env", override=False)
    model = os.environ.get("PODCAST_SCRIPT_MODEL", editorial.DEFAULT_MODEL).strip()
    fixture = json.loads((ROOT / "tests/fixtures/podcast_conversation.json").read_text(encoding="utf-8"))
    args.output.mkdir(parents=True)
    (args.output / "sources.json").write_text(json.dumps(fixture, ensure_ascii=False, indent=2), encoding="utf-8")
    manifest = {"notice": fixture["notice"], "model": model, "target_words": 180 + 280 * len(fixture["articles"]),
                "comparisons": {}, "calls": "3 calls per version, plus at most one existing caution repair"}
    current_prompts = editorial.PROMPTS
    with OpenAI() as client:
        try:
            for label, prompt_dir in (("before", args.before_prompts), ("after", current_prompts)):
                if args.only and label != args.only:
                    continue
                editorial.PROMPTS = prompt_dir
                snapshot = args.output / label
                snapshot.mkdir()
                hashes = {}
                for stage in ("plan", "write", "audit"):
                    data = (prompt_dir / f"{stage}.txt").read_text(encoding="utf-8")
                    (snapshot / f"{stage}.txt").write_text(data, encoding="utf-8")
                    hashes[stage] = hashlib.sha256(data.encode()).hexdigest()
                print(f"{label}: gerando com {model}, {manifest['target_words']} palavras-alvo", flush=True)
                try:
                    draft = editorial.generate_episode(fixture["articles"], fixture["context"], client,
                        log=lambda message: print(message, flush=True))
                except Exception as error:
                    # Do not log arbitrary SDK request bodies or secrets on failure.
                    manifest["comparisons"][label] = {"error_type": type(error).__name__, "prompt_hashes": hashes}
                    (args.output / "comparison.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
                    raise
                (snapshot / "draft.json").write_text(json.dumps(draft, ensure_ascii=False, indent=2), encoding="utf-8")
                (snapshot / "transcript.txt").write_text(editorial.transcript(draft), encoding="utf-8")
                manifest["comparisons"][label] = {"prompt_hashes": hashes, "usage": draft["usage"],
                    "rhythm": rhythm(draft["script"]), "audit": draft["audit"]}
                (args.output / "comparison.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8")
                print(json.dumps(manifest["comparisons"][label]["rhythm"], ensure_ascii=False), flush=True)
        finally:
            editorial.PROMPTS = current_prompts
    print(f"Comparação pronta: {args.output}", flush=True)


if __name__ == "__main__":
    main()
