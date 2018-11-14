from typing import Tuple,TYPE_CHECKING
if TYPE_CHECKING:
    from bulk_enumerator import BULK # type: ignore

def get_pure_struct(b:'BULK')->str:
    """
    Extract info about a prototype structure
    """
    output =  b.get_name()
    b.delete() # deallocate memory
    return output
