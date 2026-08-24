import importlib
import pkgutil
import aiida_chemshell

# Standard function to list all submodules
def list_submodules(package):
    for importer, modname, ispkg in pkgutil.walk_packages(
        package.__path__, package.__name__ + "."
    ):
        print(f"Module: {modname} | Is Sub-Package? {ispkg}")

list_submodules(aiida_chemshell)
