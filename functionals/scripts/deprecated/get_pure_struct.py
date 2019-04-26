from typing import Any
##################################################################
def get_pure_struct(b : Any)->str:
    """
    Extract info about a prototype structure
    """
    output =  b.get_name()
    b.delete() # deallocate memory
    return output
