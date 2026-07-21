import tempfile
import unittest
from pathlib import Path

from cdiff_clademl.bundle import ModelBundleError, sha256_file, verify_model_card_hash


class BundleIntegrityTests(unittest.TestCase):
    def test_matching_hash(self):
        with tempfile.TemporaryDirectory() as tmp:
            model = Path(tmp) / "model.joblib"
            model.write_bytes(b"trusted test artifact")
            digest = sha256_file(model)
            self.assertEqual(verify_model_card_hash(model, {"model_sha256": digest}), digest)

    def test_mismatched_hash_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            model = Path(tmp) / "model.joblib"
            model.write_bytes(b"changed")
            with self.assertRaisesRegex(ModelBundleError, "integrity check failed"):
                verify_model_card_hash(model, {"model_sha256": "0" * 64})


if __name__ == "__main__":
    unittest.main()
