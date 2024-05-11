"""
Translation package, will look for translation file, if exists, or creates a new one.
"""

import json
import os
import sys
from contextlib import contextmanager
from dataclasses import field
from functools import lru_cache
from pathlib import Path
from typing import Literal

type WordFormName = Literal["norm", "short", "sub", "pref"]


@lru_cache
def _deepget(path: tuple[str]):
    dictionary = get_dict()
    for key in path:
        dictionary = dictionary[key]
    return dictionary


def get_os_lang():
    os_lang = os.getenv("LANG", "en_US.UTF-8")
    lang = os_lang.split(".")[0]
    return lang


def _load_dictionary(dict_path: os.PathLike):
    with dict_path.open() as f:
        return json.load(f)


def _lang_path():
    return Path(__file__).parent / "translations"


def set_lang(lang: str):
    global __DICTIONARY
    global __LANG

    dict_path = _lang_path() / (lang + ".json")
    if not dict_path.is_file():
        print(
            f"Dict file {dict_path} is not found, switching back to en_US",
            file=sys.stderr,
        )
        return set_lang("en_US")

    _deepget.cache_clear()
    __DICTIONARY = _load_dictionary(dict_path)
    __LANG = lang


@contextmanager
def override_lang(lang: str):
    prev = get_lang()
    set_lang(lang)
    yield
    set_lang(prev)


def get_lang() -> str:
    global __LANG
    return __LANG


def available_languages() -> list[str]:
    path = _lang_path()

    return [p.stem for p in path.glob("*.json")]


def get_dict():
    global __DICTIONARY
    return __DICTIONARY


class WordForm:
    norm: str = field(compare=True)
    short: str = field(default=None, compare=False)
    sub: str = field(default=None, compare=False)

    def __init__(self, trans_path: str):
        self._trans_path = tuple(trans_path.split("."))

    @property
    def forms(self):
        forms = _deepget(self._trans_path)
        if isinstance(forms, str):
            forms = {"norm": forms}
        return forms

    def has_form(self, form: WordFormName):
        return form in self.forms

    def get(self, form: WordFormName):
        if form not in self.forms:
            return self.forms.get("norm")

        return self.forms.get(form)

    @property
    def norm(self):
        return self.get("norm")

    @property
    def short(self):
        return self.get("short")

    @property
    def sub(self):
        return self.get("sub")

    @property
    def pref(self):
        return self.get("pref")

    def __eq__(self, other):
        return isinstance(other, WordForm) and self.norm == other.norm

    def __lt__(self, other):
        return isinstance(other, WordForm) and self.norm < other.norm

    def __str__(self):
        return self.norm

    def __hash__(self):
        return hash(self.norm)

    __repr__ = __str__


def translate(key: str):
    """Translate a key to a WordForm object using the current language dictionary from buftinom/translations/"""
    return WordForm(trans_path=key)


__DICTIONARY = {}
__LANG = None
set_lang("en_US")


__all__ = [
    "translate",
    "set_lang",
    "get_lang",
    "get_dict",
    "override_lang",
    "available_languages",
    "WordForm",
    "FormParams",
]
