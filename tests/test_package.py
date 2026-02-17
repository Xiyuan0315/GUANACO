import guanaco
from guanaco import cli


def test_version_attribute_exists():
    assert isinstance(guanaco.__version__, str)
    assert guanaco.__version__


def test_cli_entrypoint_is_callable():
    assert callable(cli.main)
