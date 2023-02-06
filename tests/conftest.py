import json
import os

import pytest

server_name = "localhost"
server_port = 8123

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope="module")
def server():
    return f"http://{server_name}:{server_port}/"


@pytest.fixture(scope="module")
def off_target_body_1():
    return load_json(f"{BASE_DIR}/tests/templates/off_target_body_1.json")


@pytest.fixture(scope="module")
def off_target_body_2():
    return load_json(f"{BASE_DIR}/tests/templates/off_target_body_2.json")


@pytest.fixture(scope="module")
def off_target_body_3():
    return load_json(f"{BASE_DIR}/tests/templates/off_target_body_3.json")


@pytest.fixture(scope="module")
def off_target_body_4():
    return load_json(f"{BASE_DIR}/tests/templates/off_target_body_4.json")


@pytest.fixture(scope="module")
def off_target_body_5():
    return load_json(f"{BASE_DIR}/tests/templates/off_target_body_5.json")


@pytest.fixture(scope="module")
def off_target_body_6():
    return load_json(f"{BASE_DIR}/tests/templates/off_target_body_6.json")


@pytest.fixture(scope="module")
def off_target_body_7():
    return load_json(f"{BASE_DIR}/tests/templates/off_target_body_7.json")


@pytest.fixture(scope="module")
def off_target_body_8():
    return load_json(f"{BASE_DIR}/tests/templates/off_target_body_8.json")


@pytest.fixture(scope="module")
def off_target_body_9():
    return load_json(f"{BASE_DIR}/tests/templates/off_target_body_9.json")


@pytest.fixture(scope="module")
def on_target_body_1():
    return load_json(f"{BASE_DIR}/tests/templates/on_target_body_1.json")


@pytest.fixture(scope="module")
def on_target_body_2():
    return load_json(f"{BASE_DIR}/tests/templates/on_target_body_2.json")


@pytest.fixture(scope="module")
def flashfry_input_1():
    return load_json(f"{BASE_DIR}/tests/templates/flashfry_input_1.json")


def load_json(file_path):
    with open(file_path) as f:
        return json.load(f)
