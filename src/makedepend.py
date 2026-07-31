"""
Generates makefile directives for fortran object files with correct dependencies
based on `use` and `module` directives inside the source files
"""

import os
import sys

def make_dependencies(sources, output_file):
    module_files = {}
    file_uses = {}
    for source_file in sources:
        extension = os.path.splitext(source_file)
        file_uses[source_file] = []
        obj_file = extension[0] + '.o'
        with open(source_file) as source:
            lines = source.readlines()
            for line in lines:
                components = line.split()
                if len(components) >= 2:
                    if components[0].lower() == 'use':
                        file_uses[source_file].append(components[1])
                    if components[0].lower() == 'module':
                        module = components[1].lower()
                        if module not in module_files:
                            module_files[module] = []
                        module_files[module].append(obj_file)
    rules = []
    for module, files in module_files.items():
        module_file = module + ".mod"
        rules.append(f"{module_file}: {' '.join(files)}")

    for source_file, module_list in file_uses.items():
        extension = os.path.splitext(source_file)
        obj_file = extension[0] + '.o'
        dep_file = extension[0] + '.d'
        mod_deps = []
        for module_dep in file_uses[source_file]:
            if module_dep in module_files:
                mod_deps.append(f"{module_dep}.mod")
        rules.append(f"{obj_file}: {source_file} {' '.join(mod_deps)}")
    with open(output_file, "w") as dep:
        dep.write("\n\n".join(rules) + "\n")



if __name__ == "__main__":
    make_dependencies(sys.argv[1:], "dependencies.d")
