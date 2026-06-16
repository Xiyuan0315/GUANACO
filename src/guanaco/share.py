"""Optional public sharing for GUANACO.

Two independent pieces, both opt-in via ``settings.share``:

1. ``add_basic_auth`` — gates every HTTP request behind HTTP Basic Auth at the
   Flask layer. cloudflared quick tunnels expose the port with no access control
   of their own, so the password lives in the app, not the tunnel. This keeps the
   protection identical regardless of which tunnel tool is used.
2. ``start_cloudflared_tunnel`` — opens a cloudflared quick tunnel to the local
   port and returns the public ``https://*.trycloudflare.com`` URL.

Both degrade gracefully: a missing ``pycloudflared`` install disables tunnelling
with a hint instead of crashing the server.
"""

from __future__ import annotations

import atexit
import secrets

# Headers stamped by Cloudflare's edge on every request that arrives through a
# cloudflared tunnel. Their presence is how we tell a public visitor apart from a
# direct localhost/LAN hit (e.g. the desktop app's own webview), which carries
# none of them. A public visitor can only reach the app via the tunnel, so they
# always carry these — the gate cannot be bypassed from the internet.
_CLOUDFLARE_HEADERS = ("Cf-Ray", "Cf-Connecting-Ip")


def generate_password(n_bytes: int = 9) -> str:
    """Return a short, URL-safe random password for ad-hoc sharing."""
    return secrets.token_urlsafe(n_bytes)


def add_basic_auth(server, username: str, password: str, tunneled_only: bool = True) -> None:
    """Require HTTP Basic Auth for requests served by ``server``.

    ``server`` is the underlying Flask app (``dash.Dash().server``). Credentials
    are compared in constant time so the gate does not leak length/timing.

    With ``tunneled_only`` (the default) only requests forwarded through the
    cloudflared tunnel are challenged; direct localhost/LAN requests pass through
    untouched, so the desktop app's local webview never sees a password prompt.
    """
    from flask import Response, request

    @server.before_request
    def _require_basic_auth():  # pragma: no cover - exercised via Flask test client
        if tunneled_only and not any(h in request.headers for h in _CLOUDFLARE_HEADERS):
            return None  # direct (non-tunnel) request — not publicly shared

        auth = request.authorization
        if (
            auth is not None
            and secrets.compare_digest(auth.username or "", username)
            and secrets.compare_digest(auth.password or "", password)
        ):
            return None
        return Response(
            "GUANACO: authentication required.",
            401,
            {"WWW-Authenticate": 'Basic realm="GUANACO"'},
        )


def enable_compression(server) -> bool:
    """Enable HTTP response compression (brotli/gzip) on the Flask ``server``.

    Plotly figure payloads are JSON and compress ~5-10x, which is the difference
    between snappy and sluggish over a tunnel. Only called from ``enable_sharing``,
    so local-only runs never pay the compression CPU cost. Returns whether it was
    enabled.
    """
    try:
        from flask_compress import Compress
    except ImportError:
        print(
            "⚠️  Sharing without 'flask-compress'; responses are sent uncompressed.\n"
            "    Install it for much faster public loading:  pip install flask-compress"
        )
        return False
    Compress(server)
    return True


def enable_sharing(server, port: int, username: str, password: str) -> tuple[str | None, str, str]:
    """Gate ``server`` behind Basic Auth and open a public tunnel to ``port``.

    A blank ``password`` is replaced with a generated one. Returns
    ``(public_url, username, password)`` so the caller can show the credentials.
    """
    enable_compression(server)
    password = password or generate_password()
    add_basic_auth(server, username, password)
    public_url = start_cloudflared_tunnel(port)
    return public_url, username, password


def print_share_banner(public_url: str | None, username: str, password: str) -> None:
    """Print the public URL and login credentials (or a warning if no tunnel)."""
    print("─" * 60)
    if public_url:
        print(f"🌍 Public URL: {public_url}")
    else:
        print("🌍 Public sharing requested but no tunnel is active (see warning above).")
    print(f"🔑 Login — username: {username}   password: {password}")
    print("⚠️  Anyone with the URL AND these credentials can view your data.")
    print("─" * 60)


def start_cloudflared_tunnel(port: int) -> str | None:
    """Open a cloudflared quick tunnel to ``localhost:port``.

    Returns the public URL, or ``None`` if cloudflared is unavailable or the
    tunnel could not be started. The tunnel process is terminated at interpreter
    exit.
    """
    try:
        from pycloudflared import try_cloudflare
    except ImportError:
        print(
            "⚠️  Sharing is enabled but 'pycloudflared' is not installed.\n"
            "    Install it with:  pip install pycloudflared"
        )
        return None

    try:
        urls = try_cloudflare(port=port)
    except Exception as exc:  # binary download / network / spawn failure
        print(f"⚠️  Could not start cloudflared tunnel: {exc}")
        return None

    atexit.register(_terminate_tunnel, port)
    return getattr(urls, "tunnel", None)


def _terminate_tunnel(port: int) -> None:
    try:
        from pycloudflared import try_cloudflare

        try_cloudflare.terminate(port)
    except Exception:
        pass
