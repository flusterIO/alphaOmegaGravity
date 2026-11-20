from IPython.display import Markdown, display_markdown, display, Javascript
from IPython.core.interactiveshell import InteractiveShell
from pathlib import Path
import os
import shutil
from pyperclip import copy
import sympy as smp


def enablePrintAll():
    InteractiveShell.ast_node_interactivity = "all"  # pyright: ignore


def enablePrintLast():
    InteractiveShell.ast_node_interactivity = "last_expr_or_assign"  # pyright:ignore


def printMarkdown(content: str):
    display_markdown(Markdown(content))


def renderMathJax():
    display(
        Javascript(
            """
    let timer = null;
    (async () => {
        if (window.MathJax) {
            timer = setTimeout(() => {
            if (timer) {
                clearTimeout(timer)
            }
           await window.MathJax.typesetPromise();
                }, 50)
        }
    })();
    """
        )
    )


def clearMathImages():
    dr = Path("/Users/bigsexy/.jupyter/extensions/mathcat/mathcat/tempdir/")
    shutil.rmtree(dr)
    os.mkdir(dr)


def copyTex(content):
    copy(smp.latex(content))
    print("Copied latex successfully.")
