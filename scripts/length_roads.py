"""
@author: Elco Koks
@date: Jan, 2018
"""

import os
from functions import road_length


if __name__ == "__main__":

    base_path =   os.path.join(os.path.dirname(__file__),'..')

    road_length(base_path)