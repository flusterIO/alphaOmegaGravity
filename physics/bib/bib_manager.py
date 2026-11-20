from dataclasses import dataclass, field
import bibtexparser as bp
from bibtexparser import library
from bibtexparser.middlewares import (
    MergeCoAuthors,
    MergeNameParts,
    SplitNameParts,
    SeparateCoAuthors,
)


@dataclass
class Bib:
    file_path: str
    db: library.Library

    def __init__(
        self, file_path: str = "/Users/bigsexy/Desktop/notes/content/bib/citations.bib"
    ) -> None:
        self.file_path = file_path
        bib_database = bp.parse_file(
            file_path,
            append_middleware=[SeparateCoAuthors(), SplitNameParts()],
        )
        self.db = bib_database

    def save(self, output_path: str | None = None):
        if output_path is None:
            output_path = self.file_path
        bp.write_file(
            output_path,
            self.db,
            prepend_middleware=[MergeNameParts(), MergeCoAuthors()],
        )
