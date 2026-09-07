import ast
import io
import unittest
from contextlib import redirect_stdout
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parent


def load_functions():
    source = (ROOT / "boletim_service.py").read_text(encoding="utf-8")
    tree = ast.parse(source)
    names = {"_ler_registros_pubmed", "_ler_registros_pubmed_resiliente", "buscar_info_estruturada"}
    nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    namespace = {"time": SimpleNamespace(sleep=lambda _: None), "Entrez": None}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), "boletim_service.py", "exec"), namespace)
    return namespace


class FakeHandle:
    def close(self):
        pass


class FakeEntrez:
    def __init__(self, failures=1):
        self.failures = failures
        self.read_calls = 0
        self.fetch_calls = []

    def efetch(self, **kwargs):
        self.fetch_calls.append(kwargs["id"])
        return FakeHandle()

    def read(self, _handle):
        self.read_calls += 1
        if self.read_calls <= self.failures:
            raise RuntimeError("IncompleteRead")
        return {"PubmedArticle": []}


class PubmedFetchResilienceTests(unittest.TestCase):
    def test_retries_errors_raised_during_entrez_read(self):
        namespace = load_functions()
        fake = FakeEntrez(failures=1)
        namespace["Entrez"] = fake
        result = namespace["buscar_info_estruturada"](["1", "2"])
        self.assertEqual(result, [])
        self.assertEqual(fake.read_calls, 2)

    def test_split_failure_and_continue_without_fatal_error(self):
        namespace = load_functions()
        fake = FakeEntrez(failures=100)
        namespace["Entrez"] = fake
        with redirect_stdout(io.StringIO()):
            result = namespace["buscar_info_estruturada"](["1", "2", "3"])
        self.assertEqual(result, [])
        self.assertGreaterEqual(len(fake.fetch_calls), 3)


if __name__ == "__main__":
    unittest.main()
