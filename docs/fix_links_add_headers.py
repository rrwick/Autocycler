#!/usr/bin/env python3

import os
import re
import sys


def main():
    filename = sys.argv[1]
    lines = []
    title = os.path.basename(filename).replace('.md', '').replace('-', ' ')

    # Add a header to most pages:
    if title != 'Home' and title != 'Illustrated pipeline overview':
        lines.append(f'# {title}')
        lines.append('')

    # Load all lines and fix links:
    with open(filename, 'rt') as f:
        for line in f:
            line = line.rstrip()
            lines.append(fix_links(line))

    # Overwrite input file with new file:
    with open(filename, 'wt') as f:
        for line in lines:
            f.write(line)
            f.write('\n')


def fix_links(line):
    pattern = r'\[\[([^\|\]]+)(?:\|([^\]]+))?\]\]'
    
    def replace_link(match):
        target = match.group(2) if match.group(2) else match.group(1)  # The target (after |, or the main if no |)
        text = match.group(1) if match.group(2) else match.group(1)  # The display text (before |, or the main if no |)
        url_target = target.replace(' ', '-') + '.html'
        return f'[{text}]({url_target})'

    return re.sub(pattern, replace_link, line)


if __name__ == '__main__':
    main()


def test_fix_links():
    assert fix_links('[[Generating input assemblies]]') == '[Generating input assemblies](Generating-input-assemblies.html)'
    assert fix_links('[[generate input assemblies for Autocycler|Generating input assemblies]]') == '[generate input assemblies for Autocycler](Generating-input-assemblies.html)'
    assert fix_links('Lorem [[ipsum]] dolor [[sit]] amet') == 'Lorem [ipsum](ipsum.html) dolor [sit](sit.html) amet'
    assert fix_links('Lorem [[ipsum|i]] dolor [[sit|s]] amet') == 'Lorem [ipsum](i.html) dolor [sit](s.html) amet'

