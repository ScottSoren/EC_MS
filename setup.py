import codecs
import os
import re
from setuptools import setup, find_packages

NAME = "EC_MS"
KEYWORDS = [
    "electrochemistry",
    "mass spectrometry",
]

CLASSIFIERS = [
    "Intended Audience :: Scientists",
    "Natural Language :: English",
    "License :: MIT License",
]

PACKAGES = find_packages(where="src")

# print(PACKAGES)  # test

HERE = os.path.abspath(os.path.dirname(__file__))
META_PATH = os.path.join("src", "EC_MS", "__init__.py")


def read(*parts):
    """
    Build an absolute path from *parts* and return the contents of the
    resulting file. Assume UTF-8 encoding.
    """
    with codecs.open(os.path.join(HERE, *parts), "rb", "utf-8") as f:
        return f.read()


META_FILE = read(META_PATH)
INSTALL_REQUIRES = []


def find_meta(meta):
    """
    Extract __*meta*__ from META_FILE
    """
    meta_match = re.search(
        r"^__{meta}__ = ['\"]([^'\"]*)['\"]".format(meta=meta), META_FILE, re.M
    )
    if meta_match:
        print("Able to find __{meta}__ string".format(meta=meta))
        return meta_match.group(1)
    print("Unable to find __{meta}__ string".format(meta=meta))
    raise RuntimeError("Unable to find __{meta}__ string".format(meta=meta))


# version = find_meta("version"),
# url = find_meta("url"),
# author = find_meta("author"),
# print("{}\n{}\n{}".format(version, url, author))


if __name__ == "__main__":
    setup(
        name=NAME,
        description=find_meta("description"),
        license=find_meta("license"),
        version=find_meta("version"),
        url=find_meta("url"),
        author=find_meta("author"),
        author_email=find_meta("email"),
        maintainer=find_meta("author"),
        maintainer_email=find_meta("email"),
        keywords=KEYWORDS,
        packages=PACKAGES,
        package_dir={"": "src"},
        zip_safe=False,
        include_package_data=True,
        classifiers=CLASSIFIERS,
        install_requires=INSTALL_REQUIRES,
        options={"bdist_wheel": {"universal": "1"}},
    )
