import re
import sys
from pathlib import Path


TABLE_NAMES = ("Table1", "Table2", "Table3")


def replace_table_file(text, table_name, table_path):
    block_pattern = re.compile(
        rf"(?P<header>^[ \t]*{re.escape(table_name)}[ \t]*\n"
        rf"[ \t]*\{{[ \t]*\n)"
        rf"(?P<body>.*?)"
        rf"(?P<footer>^[ \t]*\}}[ \t]*$)",
        flags=re.MULTILINE | re.DOTALL,
    )
    block_matches = list(block_pattern.finditer(text))
    if len(block_matches) != 1:
        raise ValueError(
            f"expected exactly one {table_name} dictionary, "
            f"found {len(block_matches)}"
        )

    match = block_matches[0]
    file_pattern = re.compile(
        r'^(?P<prefix>[ \t]*file[ \t]+)(?:"[^"\n]*"|[^;\n]*)(?P<suffix>;.*)$',
        flags=re.MULTILINE,
    )
    body, replacement_count = file_pattern.subn(
        lambda file_match: (
            f'{file_match.group("prefix")}"{table_path}"'
            f'{file_match.group("suffix")}'
        ),
        match.group("body"),
        count=1,
    )
    if replacement_count != 1:
        raise ValueError(f"expected exactly one file entry for {table_name}")

    replacement = match.group("header") + body + match.group("footer")
    return text[: match.start()] + replacement + text[match.end() :]


def configure_table_files(properties_path, table_paths):
    if len(table_paths) != len(TABLE_NAMES):
        raise ValueError(f"expected {len(TABLE_NAMES)} table paths")

    properties_path = Path(properties_path)
    text = properties_path.read_text(encoding="utf-8")
    for table_name, table_path in zip(TABLE_NAMES, table_paths):
        text = replace_table_file(text, table_name, table_path)
    properties_path.write_text(text, encoding="utf-8")


if __name__ == "__main__":
    if len(sys.argv) != 5:
        raise SystemExit(
            "usage: configureTableFiles.py COMBUSTION_PROPERTIES "
            "TABLE1_PATH TABLE2_PATH TABLE3_PATH"
        )
    configure_table_files(Path(sys.argv[1]), sys.argv[2:])
