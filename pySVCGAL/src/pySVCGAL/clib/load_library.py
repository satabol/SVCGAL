import ctypes
import os 
from sys import platform 

#print("==== pySVCGAL ==========================================")

def load_library():
    #print("==== load_library pySVCGAL ==========================================")
    if platform == 'linux' or platform == 'linux2':        
        SVCGAL_clib = ctypes.CDLL(os.path.join(os.path.dirname(__file__),'libctSVCGAL.so'))
        pass
    elif platform == 'darwin':
        # OSX
        SVCGAL_clib = ctypes.CDLL(os.path.join(os.path.dirname(__file__),'libctSVCGAL.dylib'))
        pass
    elif platform == 'win32':        
        here = os.path.dirname(__file__).replace('\\','/') 
        SVCGAL_clib = ctypes.CDLL(os.path.join(here,"ctSVCGAL.dll"))
    return SVCGAL_clib

if __name__ == "__main__":
    SVCGAL_clib = load_library()