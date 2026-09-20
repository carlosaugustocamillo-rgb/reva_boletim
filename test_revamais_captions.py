import ast
import re
import unittest
from pathlib import Path


def load_caption_helpers():
    source = Path("revamais_service.py").read_text(encoding="utf-8")
    module = ast.parse(source)
    needed = {"resumir_texto_para_legenda", "legenda_repete_texto", "normalizar_legenda_visual"}
    functions = [
        node for node in module.body
        if isinstance(node, ast.FunctionDef) and node.name in needed
    ]
    namespace = {"re": re}
    exec(compile(ast.Module(body=functions, type_ignores=[]), "caption_helpers", "exec"), namespace)
    return namespace


class RevaMaisCaptionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.helpers = load_caption_helpers()

    def test_rejects_copy_of_next_paragraph(self):
        section = (
            "A caminhada supervisionada melhora a tolerância ao esforço em pessoas com limitação funcional. "
            "Os resultados devem ser interpretados de acordo com a população estudada."
        )
        caption = "A caminhada supervisionada melhora a tolerância ao esforço em pessoas com limitação funcional."
        self.assertTrue(self.helpers["legenda_repete_texto"](caption, section))
        self.assertEqual(self.helpers["normalizar_legenda_visual"](caption, section), "")

    def test_keeps_visual_interpretation_that_adds_context_without_copying(self):
        section = (
            "A caminhada supervisionada melhora a tolerância ao esforço em pessoas com limitação funcional."
        )
        caption = "A progressão gradual representa adaptação funcional acompanhada pelo fisioterapeuta."
        self.assertFalse(self.helpers["legenda_repete_texto"](caption, section))
        self.assertEqual(
            self.helpers["normalizar_legenda_visual"](caption, section),
            caption,
        )


if __name__ == "__main__":
    unittest.main()
