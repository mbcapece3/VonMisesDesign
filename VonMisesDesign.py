from VMPlot import VMPlot
from VonMises import VonMises

titlestring = """
 __     __          __  __ _               ____            _             
 \ \   / /__  _ __ |  \/  (_)___  ___  ___|  _ \  ___  ___(_) __ _ _ __  
  \ \ / / _ \| '_ \| |\/| | / __|/ _ \/ __| | | |/ _ \/ __| |/ _` | '_ \ 
   \ V / (_) | | | | |  | | \__ \  __/\__ \ |_| |  __/\__ \ | (_| | | | |
    \_/ \___/|_| |_|_|  |_|_|___/\___||___/____/ \___||___/_|\__, |_| |_|
                                                             |___/       
"""

print(titlestring)
print("Developed by Matthew Capece \nVonMisesDesign v1.0 \nMarch 23, 2024 \n\n\n")

center = 0+0j
singularities = [1+0j, -1+0j]

v_inf = 1
alpha_rad = 0

vm1 = VonMises(center, v_inf, alpha_rad, singularities)

vmPlot = VMPlot(vm1)