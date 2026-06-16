import base64

from flask import Flask

from guanaco import share


def _client_with_auth(username="guanaco", password="secret"):
    app = Flask(__name__)

    @app.route("/")
    def index():
        return "protected"

    share.add_basic_auth(app, username, password)
    return app.test_client()


# A cloudflared-forwarded (public) request carries a Cf-Ray header; a direct
# localhost/LAN request (e.g. the desktop webview) does not.
TUNNEL_HEADER = {"Cf-Ray": "abc123-SJC"}


def _auth_header(username, password, *, tunneled=True):
    token = base64.b64encode(f"{username}:{password}".encode()).decode()
    headers = {"Authorization": f"Basic {token}"}
    if tunneled:
        headers.update(TUNNEL_HEADER)
    return headers


def test_direct_request_passes_without_credentials():
    # No Cloudflare header => local/direct request => not challenged.
    resp = _client_with_auth().get("/")
    assert resp.status_code == 200
    assert resp.get_data(as_text=True) == "protected"


def test_tunneled_request_without_credentials_is_rejected():
    resp = _client_with_auth().get("/", headers=TUNNEL_HEADER)
    assert resp.status_code == 401
    assert "Basic" in resp.headers["WWW-Authenticate"]


def test_tunneled_request_with_wrong_password_is_rejected():
    resp = _client_with_auth().get("/", headers=_auth_header("guanaco", "wrong"))
    assert resp.status_code == 401


def test_tunneled_request_with_correct_credentials_passes():
    resp = _client_with_auth().get("/", headers=_auth_header("guanaco", "secret"))
    assert resp.status_code == 200
    assert resp.get_data(as_text=True) == "protected"


def test_generated_password_is_nonempty_and_unique():
    p1 = share.generate_password()
    p2 = share.generate_password()
    assert p1 and p2
    assert p1 != p2


def test_enable_sharing_generates_blank_password_and_gates(monkeypatch):
    # Don't actually open a tunnel during the test.
    monkeypatch.setattr(share, "start_cloudflared_tunnel", lambda port: None)

    app = Flask(__name__)

    @app.route("/")
    def index():
        return "protected"

    url, username, password = share.enable_sharing(app, 4399, "guanaco", "")
    assert url is None
    assert username == "guanaco"
    assert password  # a random password was generated for the blank input

    client = app.test_client()
    assert client.get("/", headers=TUNNEL_HEADER).status_code == 401
    assert client.get("/", headers=_auth_header("guanaco", password)).status_code == 200


def test_enable_compression_degrades_without_flask_compress(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "flask_compress":
            raise ImportError("not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    assert share.enable_compression(Flask(__name__)) is False


def test_tunnel_returns_none_without_pycloudflared(monkeypatch):
    # Simulate pycloudflared not being installed: the import inside the function fails.
    import builtins

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "pycloudflared":
            raise ImportError("not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    assert share.start_cloudflared_tunnel(4399) is None
