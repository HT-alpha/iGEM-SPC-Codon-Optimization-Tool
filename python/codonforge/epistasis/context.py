"""Install `epistasis.context_naturalness` from the verified Potts artifact.

The couplings stay in the process cache. A tier request carries the digest.
A missing artifact is a refusal, not a score of zero. The dwell head is a
different term and is not loaded here.
"""

from __future__ import annotations

import json
from typing import Any

from codonforge import _core
from codonforge.epistasis.artifact import artifact_path, load_calibration, verify_artifact
from codonforge.modelreg.store import ArtifactStore

TERM_ID = "epistasis.context_naturalness"


class ContextError(Exception):
    """The context term cannot be scored. Carries a `CF-*` code."""

    def __init__(self, code: str, detail: str) -> None:
        self.code = code
        self.detail = detail
        super().__init__(f"{code}: {detail}")


def prepare(store: ArtifactStore, host_id: str) -> dict[str, Any]:
    """Publish the public E. coli artifact and install its selected model.

    The artifact's host is MG1655. A job whose host id starts with
    `escherichia_coli:` may pin it. No distance penalty is applied, because
    none has been measured. Any other host is refused.
    """
    if not str(host_id).startswith("escherichia_coli:"):
        raise ContextError(
            "CF-LAND-001",
            "epistasis.context_naturalness has a verified artifact for "
            "escherichia_coli:511145/MG1655 only; "
            f"{host_id} is not that host, and no other Potts artifact is installed",
        )
    path = artifact_path()
    if not path.is_file():
        raise ContextError(
            "CF-LAND-001",
            "epistasis.context_naturalness was requested but the Potts artifact "
            f"is absent at {path}. The term is not scored as zero",
        )
    document = json.loads(path.read_text(encoding="utf-8"))
    return install_document(store, document)


def install_pinned(store: ArtifactStore, digest: str, artifact_id: str) -> dict[str, Any]:
    """Install the model a sealed lockfile already names.

    The store is the first place looked. The public file is used only when
    its digest is the one the lockfile recorded.
    """
    manifest = store.get(digest)
    if manifest is None:
        path = artifact_path()
        if path.is_file():
            document = json.loads(path.read_text(encoding="utf-8"))
            verify_artifact(document)
            if document.get("digest") == digest:
                manifest = store.publish(document)
        if manifest is None:
            raise ContextError(
                "CF-RESOLVE-006",
                f"pinned artifact digest {digest} is not in the store and the "
                "public Potts file does not match it; replay cannot proceed",
            )
    if manifest.artifact_id != artifact_id:
        raise ContextError(
            "CF-RESOLVE-006",
            f"pin names {artifact_id} but digest {digest} holds {manifest.artifact_id}",
        )
    return _install(manifest.document)


def install_document(store: ArtifactStore, document: dict[str, Any]) -> dict[str, Any]:
    verify_artifact(document)
    manifest = store.publish(document)
    if manifest.digest != document["digest"]:
        raise ContextError(
            "CF-RESOLVE-006",
            f"the Potts artifact records digest {document['digest']} but publishing "
            f"it produced {manifest.digest}",
        )
    return _install(manifest.document)


def _install(document: dict[str, Any]) -> dict[str, Any]:
    parameters = document["parameters"]
    window = int(parameters["selected_window"])
    chosen = next(
        (row for row in parameters["models"] if int(row["window"]) == window),
        None,
    )
    if chosen is None:
        raise ContextError(
            "CF-LAND-001",
            f"the Potts artifact selects window {window} but does not contain that model",
        )
    payload = {
        "digest": document["digest"],
        "artifact_id": document["artifact_id"],
        "window": window,
        "h": chosen["h"],
        "couplings": chosen["couplings"],
    }
    json.loads(_core.install_context_model(json.dumps(payload)))
    # load_calibration is what the page is allowed to say. Touching it here
    # keeps a refused artifact from being described as fitted.
    calibration = load_calibration()
    if calibration.get("state") != "fitted":
        raise ContextError(
            "CF-LAND-001",
            "epistasis.context_naturalness artifact did not verify as fitted",
        )
    if calibration.get("digest") != document["digest"]:
        raise ContextError(
            "CF-RESOLVE-006",
            "the installed Potts digest does not match the calibration document",
        )
    return {
        "digest": document["digest"],
        "artifact_id": document["artifact_id"],
        "window": window,
    }
