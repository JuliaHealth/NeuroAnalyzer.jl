#!/usr/bin/env python3
"""Publish a new version of NeuroAnalyzer.jl on zenodo.org.

A new version of the latest Zenodo record is created; its files are replaced with the release
ZIP archive, its version and publication date are updated (all other metadata are copied from
the previous version) and it is published.

usage: publish_zenodo.py [--draft] [--check] <version> <zip file>
  --draft  do not publish: leave the new version as a draft for review on zenodo.org
  --check  only check the token (no version or ZIP file needed)

environment:
  NA_ZENODO_TOKEN  zenodo.org personal access token with the deposit:write and deposit:actions
                   scopes

Re-running is safe: a version that is already published is skipped and an unpublished draft
left by a previous run is reused.
"""

import datetime
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

API = "https://zenodo.org/api"
# concept record: all NeuroAnalyzer versions (DOI 10.5281/zenodo.7372648)
CONCEPT_ID = "7372648"


def fail(msg):
    print(f"Error: {msg}", file=sys.stderr)
    sys.exit(1)


def request(method, url, token=None, data=None, json_data=None, headers=None):
    hdrs = dict(headers or {})
    if token:
        hdrs["Authorization"] = f"Bearer {token}"
    if json_data is not None:
        data = json.dumps(json_data).encode()
        hdrs["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, method=method, headers=hdrs)
    try:
        with urllib.request.urlopen(req, timeout=600) as resp:
            body = resp.read()
    except urllib.error.HTTPError as e:
        fail(f"{method} {url}: HTTP {e.code}: {e.read().decode(errors='replace')[:500]}")
    except urllib.error.URLError as e:
        fail(f"{method} {url}: {e.reason}")
    return json.loads(body) if body else None


def main(args):
    draft_only = "--draft" in args
    check_only = "--check" in args
    pos = [a for a in args if not a.startswith("--")]
    unknown = [a for a in args if a.startswith("--") and a not in ("--draft", "--check")]
    if unknown or (not check_only and len(pos) != 2):
        fail("usage: publish_zenodo.py [--draft] [--check] <version> <zip file>")

    token = os.environ.get("NA_ZENODO_TOKEN")
    if not token:
        fail("NA_ZENODO_TOKEN not set")

    # the token must be valid (listing own depositions needs the deposit scopes)
    request("GET", f"{API}/deposit/depositions?size=1", token)
    if check_only:
        print("zenodo.org token: OK")
        return

    version, zip_file = pos
    if not os.path.isfile(zip_file):
        fail(f"file not found: {zip_file}")

    # latest published version
    latest = request("GET", f"{API}/records/{CONCEPT_ID}/versions/latest")
    if latest["metadata"].get("version") == version:
        print(f"Version {version} is already published: {latest['links']['html']}")
        return

    # new version (an existing unpublished draft is returned instead of a new one)
    dep = request("POST", f"{API}/deposit/depositions/{latest['id']}/actions/newversion", token)
    draft = request("GET", dep["links"]["latest_draft"], token)
    draft_id = draft["id"]
    print(f"Draft: https://zenodo.org/deposit/{draft_id}")

    # replace the files copied from the previous version
    for f in draft.get("files", []):
        request("DELETE", f"{API}/deposit/depositions/{draft_id}/files/{f['id']}", token)
    name = os.path.basename(zip_file)
    print(f"Uploading {name}")
    with open(zip_file, "rb") as fh:
        request(
            "PUT",
            f"{draft['links']['bucket']}/{urllib.parse.quote(name)}",
            token,
            data=fh,
            headers={
                "Content-Type": "application/octet-stream",
                "Content-Length": str(os.path.getsize(zip_file)),
            },
        )

    # metadata of the previous version with the new version number and date;
    # the DOI of the draft is reserved by zenodo
    meta = draft["metadata"]
    meta.pop("doi", None)
    meta["version"] = version
    meta["publication_date"] = datetime.date.today().isoformat()
    request("PUT", f"{API}/deposit/depositions/{draft_id}", token, json_data={"metadata": meta})

    if draft_only:
        print(f"Draft ready for review (not published): https://zenodo.org/deposit/{draft_id}")
        return

    pub = request("POST", f"{API}/deposit/depositions/{draft_id}/actions/publish", token)
    print(f"Published version {version}: DOI {pub.get('doi')}, {pub['links'].get('html')}")


if __name__ == "__main__":
    main(sys.argv[1:])
