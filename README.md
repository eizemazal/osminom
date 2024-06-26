# OSMINOM - Open source SMILES and IUPAC name converter #

## Synopsis ##
Цель проекта - демонстрация возможности создания конвертера названий ИЮПАК в структуры в пределах 500 строк на Python на основе контекстно свободных грамматик с одной зависимостью от lex / yacc на Python (ply).
Парсинг в этой версии написан для русского языка, его очень легко сделать многоязычным, но сейчас этого не делаю в пользу легкости чтения кода.


## Запуск ##

Нужен Python 3.10+. Настройка виртуального окружения:
```
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

```

Запуск тестов
```
python3 -m pytest tests/
```

## Примеры работы ##
- (2-(метиламино)этил)бензол → CNCCc1ccccc1
- 1-гидрокси-3-метилбензол → c1ccc(C)cc1O
- 1-метил-1-(2-гидроксифенил)-2-цианоциклобутан → C1CC(C#N)C1(C)c1cccc(c1)O


## Токены и правила ##

TBASE - тривиальное название родительского соединения. Бензол, нафталин и т.д.

ALIPHATIC - корень алифатической цепи. Мет, эт, проп, бут, ...

alicyclic - цикло ALIPHATIC >= 3. Циклопроп, циклобут, ...

SUFFIX - суффикс названия, который совместим с корнями цепей. Ан, ен, ин, аль, ановая кислота.

TODO : SUFFIX_T - терминальный суффикс, после которого нельзя размещать суффиксные формы. ановая кислота, ансульфокислота.

base (основа) - цепь с суффиксом либо тривиальное название.
Метан, бутановая кислота. Это уже годится как полное название соединения.

lform - (3-метилбутил), гидрокси, метокси, ((4-(4-аминофенил)фенил)этил)метиламино. Перед основой может быть несколько лформ с локантами или без них.

llocant - локант слева. 5-, 2-, 12-

PRFXFORM - форма префикса заместителя, в которой нельзя делать замещение.

TRADICAL - префиксы заместителей, в которых можно делать замещение, например, аллил, фенил, амино.

## ToDo - маленькие задачи

- доработать поддержку существующего кода, префиксных форм, дописать тесты
- добавить поддержку суффиксных форм (бутанол-2, бутан-2-ол-1-амин, пентановая кислота, хлорциклогексен-2-карбальдегид)
- поддержка префиксов типа 5,5-диметил, 1,2,3-тригидрокси
- поддержка непредельных суффиксов ен и ин, которые нельзя записать как прибавляемые к родительскому фрагменты (переписать прибавление как объект-оператор, добавить оператор для дегидрирования и других модификаций, например замены атома)
- 8-окса, 9-аза и подобные замены углерода на гетероатомы
- добавить многоязычность с поддержкой английского языка
- добавить поддержку стереохимии двойной связи - вывод SMILES со стереохимией, поддержку Z,E и цис-транс номенклатуры
- добавить поддержку RS-стереохимии в локанты и конвертацию ее в стереохимию SMILES
- добавить больше родительских соединений в TBASE, суффиксов и префиксов

## ToDo - большие задачи
- добавить поддержку гетероциклических, би- и полициклических конденсированных и каркасных соединений. Написать тесты, дизазабицикло[5.4.0]ундецен-7, пентацикло[4.2.0.02,5.03,8.04,7]октан (кубан), пурин и 1,3-тиазол должны проходить.
- добавить конвертацию структур в названия
