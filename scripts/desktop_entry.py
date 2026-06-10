"""Compiled entry point for the GUANACO desktop app.

Nuitka compiles this thin launcher into the standalone binary / .app bundle.
Keeping it tiny (just a call into the package) means all the real logic stays
in ``guanaco.desktop`` and is importable/testable normally.
"""

from guanaco.desktop import desktop_main

if __name__ == "__main__":
    desktop_main()
