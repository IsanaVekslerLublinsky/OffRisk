import httpx
import pytest
import json


def test_get_root(server):
    r = httpx.get(server)
    assert r.status_code == 200
    assert r.text == "Welcome to off-tov server"


@pytest.mark.parametrize("run_type, body", [("off", "off_target_body_1"),
                                            ("off", "off_target_body_2"),
                                            ("off", "off_target_body_3"),
                                            ("off", "off_target_body_4"),
                                            ("off", "off_target_body_5"),
                                            ("off", "off_target_body_6"),
                                            ("off", "off_target_body_7"),
                                            ("off", "off_target_body_8"),
                                            ("off", "off_target_body_9"),
                                            ("off", "off_target_body_10"),
                                            ("off", "off_target_body_11"),
                                            ("off", "off_target_body_12"),
                                            ("off", "off_target_body_13"),
                                            ("off", "off_target_body_14"),
                                            ("on", "on_target_body_1"),
                                            ("on", "on_target_body_2")])
def test_off_target(server, run_type, body, request):
    r = httpx.post("{}/v1/{}-target-analyze/".format(server, run_type), timeout=10000,
                   json=request.getfixturevalue(body))
    assert r.status_code == 200
    request_response = r.json()
    keys = list(request_response.keys())
    assert keys == ['request_id', 'off_targets', 'flashfry_score', 'all_result', 'time']


@pytest.mark.parametrize("dna_bulge, rna_bulge, run_type", [("0", "0", "cas-offinder-bulge"),
                                                           ("2", "4", "cas-offinder-bulge"),
                                                           ("0", "0", "cas-offinder")])
def test_cas_offinder(server, dna_bulge, rna_bulge, run_type):
    r = httpx.get("{}/v1/{}/"
                  "?genome=human&pattern=NNNNNNNNNNNNNNNNNNNNNGG%20{}%20{}&sequences=GTTTTTTGTTGACCCGGAAACGG%204"
                  .format(server, run_type, dna_bulge, rna_bulge), timeout=10000)
    assert r.status_code == 200
    request_response = json.loads(r.text)
    keys = list(request_response.keys())
    assert keys == ["cas-offinder-dataframe", "message"]


@pytest.mark.parametrize("body", ["flashfry_input_1"])
def test_flashfry(server, body, request):
    r = httpx.post("{}/v1/flashfry/".format(server), timeout=10000,
                   json=request.getfixturevalue(body))
    assert r.status_code == 200
    request_response = r.json()
    keys = list(request_response.keys())
    assert keys == ["flashfry-discover", "flashfry-score", "message"]
