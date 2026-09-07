"""Local and Plotly Cloud entry point: python examples/linked_views/app.py."""

import os

from showcase import create_app

app = create_app()
server = app.server

if __name__ == "__main__":
    app.run(host="127.0.0.1", port=int(os.environ.get("PORT", "8050")), debug=False)
