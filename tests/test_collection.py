import importlib.util
from pathlib import Path
import tempfile
import unittest
import zipfile
spec=importlib.util.spec_from_file_location('collector',Path(__file__).resolve().parents[1]/'tools/collect_production_sources.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class SourceCollection(unittest.TestCase):
    def test_only_sources_and_no_replacement(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/'production';root.mkdir()
            (root/'producer.C').write_text('void producer() {}\n')
            (root/'job.sub').write_text('queue 1\n')
            (root/'events.root').write_bytes(b'data')
            (root/'.env').write_text('local-settings')
            (root/'results').mkdir();(root/'results/old.py').write_text('x=1')
            (root/'linked.C').symlink_to(root/'producer.C')
            _,rows=module.collect(root)
            self.assertEqual({r['path'] for r in rows},{'producer.C','job.sub'})
            dest=Path(tmp)/'sources.zip';module.collect(root,dest)
            with zipfile.ZipFile(dest) as z:
                self.assertEqual(len(z.namelist()),3)
                self.assertIn('production_sources/SOURCE_INVENTORY.json',z.namelist())
            before=dest.read_bytes()
            with self.assertRaises(FileExistsError):module.collect(root,dest)
            self.assertEqual(dest.read_bytes(),before)
