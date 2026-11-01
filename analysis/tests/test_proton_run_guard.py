import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow"))
from proton_selection import guard_proton_run


class ProtonRunGuard(unittest.TestCase):
    def test_no_reuse_of_old_outputs_or_mixed_policy(self):
        with tempfile.TemporaryDirectory() as d, patch.dict(os.environ, {"YPSILON_PROTON_SELECTION": "required"}):
            root=Path(d); old=root/'old'; old.mkdir(); output=old/'result.csv'; output.write_bytes(b'old output')
            config={'run': {'dir':str(old)}}
            with self.assertRaises(ValueError): guard_proton_run(config,root)
            self.assertEqual(output.read_bytes(),b'old output')
            self.assertFalse((old/'proton_selection.json').exists())
            fresh=root/'fresh'; config={'run': {'dir':str(fresh)}}
            guard_proton_run(config,root,dry_run=True); self.assertFalse(fresh.exists())
            guard_proton_run(config,root); guard_proton_run(config,root)
            self.assertEqual(config['analysis']['proton_selection']['mode'],'required')
            with patch.dict(os.environ, {'YPSILON_PROTON_SELECTION':'inclusive'}), self.assertRaises(ValueError):
                guard_proton_run(config,root)
            source=root/'scripts/a.py';source.parent.mkdir();source.write_text('x=1')
            with self.assertRaises(ValueError): guard_proton_run(config,root)

if __name__=='__main__': unittest.main()
