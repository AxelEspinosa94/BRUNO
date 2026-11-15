import os

from app.utils.LayoutEngine import LayoutEngine

base = os.path.dirname(__file__)

lE = LayoutEngine(base=base, module="ode_ivp")

ode_ivp_section = lE.build()
