"""
Translation package, will look for translation file, if exists, or creates a new one.
"""

import json
import os
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import KeysView, Literal, TypedDict, Unpack


@dataclass(frozen=True, eq=True, unsafe_hash=True, order=True)
class WordForm:
    norm: str = field(compare=True)
    short: str = field(default=None, compare=False)

    def get(self, form: Literal["norm", "short"]):
        return getattr(self, form)

    def __str__(self):
        return self.norm

    __repr__ = __str__


class WordFromDict(TypedDict):
    norm: str
    short: str


L = KeysView[WordFromDict]


class FormParams(TypedDict):
    short: bool


_REGISTERED_TRANSLATIONS = {}


def _deepset(dictionary: dict, path: list, value):
    for key in path[:-1]:
        dictionary = dictionary.setdefault(key, {})
    dictionary[path[-1]] = value


def _deepget(dictionary: dict, path: list):
    for key in path:
        dictionary = dictionary[key]
    return dictionary


def translate(key: str, /, **forms: Unpack[FormParams]) -> WordForm:
    """Translate a key to a WordForm object using the current language dictionary from buftinom/translations/"""
    """Will be overriden by the actual translation function, this is just a stub for correct type hints"""
    pass


def _translate(dictionary, key, /, **forms: Unpack[FormParams]):
    formDict = _deepget(dictionary, key.split("."))
    if not forms and isinstance(formDict, str):
        return WordForm(norm=formDict)

    for form, required in forms.items():
        if required and not formDict.get(form):
            raise ValueError(f"Missing required form {form} for {key}")
    return WordForm(**formDict)


def _translate_tmp(word, /, **forms: Unpack[FormParams]):
    return WordForm(norm=word, short=word)


def _translate_register(dict_path: os.PathLike, key, /, **forms: Unpack[FormParams]):
    print(f"Registering translation for {key}")
    fq_path = key.split(".")
    word = fq_path[-1]

    if key not in _REGISTERED_TRANSLATIONS:
        _deepset(
            _REGISTERED_TRANSLATIONS,
            fq_path,
            {
                "norm": word,
                **{form: word for form, required in forms.items() if required},
            },
        )

        with dict_path.open("w") as f:
            json.dump(_REGISTERED_TRANSLATIONS, f, indent=4)

    return _translate_tmp(word, **forms)


def _load_dictionary(dict_path: os.PathLike):
    if dict_path.exists():
        with dict_path.open() as f:
            return json.load(f)

    return {}


def _get_lang():
    os_lang = os.getenv("LANG", "en_US.UTF-8")
    lang = os_lang.split(".")[0]
    return lang


def _resolve_translate_func(mode, lang):
    if mode:
        return _translate_register

    dict_path = Path(__file__).parent / f"translations/{lang}.json"
    dictionary = _load_dictionary(dict_path)
    if not dictionary:
        return partial(_translate_register, dict_path)
    else:
        return partial(_translate, dictionary)


translate = _resolve_translate_func(
    os.getenv("__BUILDS_TRANSLATIONS", False),
    _get_lang(),
)

__all__ = [
    "translate",
    "WordForm",
    "FormParams",
]
