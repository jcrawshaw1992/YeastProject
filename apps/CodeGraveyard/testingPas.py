import argparse
parser = argparse.ArgumentParser()
parser.add_argument('foo')
parser.add_argument('bar', nargs='?', default='eggs')
with assertRaisesRegex(ArgumentParseError) as cm:
       parser.parse_args([])
   self.assertNotIn('bar', str(cm.exception))