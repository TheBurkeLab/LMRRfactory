import LMRRfactory
from importlib import resources
import yaml, pathlib

print("Version:", LMRRfactory.__version__)
print("Installed at:", LMRRfactory.__file__)

# Adjust the path to match where the YAML *should* be inside the package
yaml_path = pathlib.Path(LMRRfactory.__file__).parent / "data" / "db.yaml"
print("YAML path:", yaml_path)
print("YAML exists:", yaml_path.exists())

if yaml_path.exists():
    print("YAML contents as seen by interpreter:")
    print(yaml_path.read_text())