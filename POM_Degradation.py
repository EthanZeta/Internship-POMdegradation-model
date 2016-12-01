
import numpy as np
import matplotlib.pyplot as plt

R = 8.314  # R constant for ideal gas, [J/(mol*K)]


# input several conditions and constant for further calculation
def partial_pressure(fraction=3, tem=283.15, pre=1.21325):
    """
    This function is used to define basic conditions for calculation
    :param fraction: float, HCl gas fraction in pipe, [ppm]
    :param tem: float, temperature in pipe, [K]
    :param pre: float, absolute pressure in the pipe, [Bar]
    :return: float, temperature [K] and pressure [Bar] values for other functions, partial pressure of water vapour [Bar] is
    computed by vapour formula in air, p_hcl [Bar] is the partial pressure value in the pipe
    """
    pwater = np.exp(77.345 + 0.0057 * tem - 7235 / tem) / (tem ** 8.2) / 100000
    # partial pressure of H2O in air by relation, [Bar]
    p_hcl = fraction * 10 ** -5 * pre
    # firstly use 3ppm concentration to do estimation [Bar]
    return tem, pre, pwater, p_hcl


# calculate the concentration of each gas at the surface and bulk
def average_concentration():
    """
    This function is used to calculate the fraction value of all the components.
    :return: float, xid is fraction difference between droplet surface and the bulk,
    xi_bar is the average fraction in the thin film
    """
    x10 = 0     # HCl gas fraction on the droplet surface, [/]
    x30 = pwater / pre      # water vapour fraction in the pipe, [/]
    x11 = phcl / pre        # HCl gas fraction in the pipe, [/]
    x31 = 12.3e-3     # water vapour fraction on the droplet surface, [/]
    if x30 < 0.02 and x31 < 0.02:
        x30 = 0.0
        x31 = 0.0       # the fraction is pretty low, thus neglect the water part
    x20 = 1 - x10 - x30     # other gas fraction on the droplet surface, [/]
    x21 = 1 - x11 - x31     # other gas fraction in the pipe, [/]
    x1d = x10 - x11     # HCl fraction difference, [/]
    x1_bar = (x10 + x11) / 2  # HCl average fraction, [/]
    x2d = x20 - x21  # Natural gas fraction difference, [/]
    x2_bar = (x20 + x21) / 2    # natural gas average fraction, [/]
    x3d = x30 - x31  # water fraction difference, [/]
    x3_bar = (x30 + x31) / 2  # water average fraction, [/]
    return x1d, x1_bar, x2d, x2_bar, x3d, x3_bar


# Diffusion coefficient of each component
def flux():
    """
    Using the thin film approximation and Maxwell-Stefan equations, we can calculate the HCl flux on the water surface
    :return: float, N1 is the HCl flux, converted into the same units for further calculation, [mol/(dm2*sec)]
    """
    delta = 0.01  # film thickness, [dm]
    c = pre * 10 ** 2 / (R * tem)  # total concentration calculated by ideal gas equation, in [mol/L]
    D12 = 0.001626528 / pre   # HCl diffusion in Air, [dm2/s] @296K
    D13 = 3e-7  # HCl gas diffusion in water, [dm2/s] @296K
    D23 = 1.5e-7  # CH4 gas diffusion in water, [dm2/s] @296K
    N1 = ((x1_bar * x2d * D23) / (x2_bar * delta * D13) - x1_bar / delta) / \
         (x2_bar / (D12 * c) + x3_bar / (D13 * c) + D23 * x1_bar / (D12 * D13 * c))
    # print 'Flux of HCl into water', abs(N1), [mol/(dm2*sec)]
    return N1


# calculation for the condensed water volume and surface as droplet attached on a flat surface
# (assume the bottom of pipe can be accounted as flat)
def droplet(r_drop=0.02):  # [dm]
    """
    The function is calculating the geometric parameters of single droplet.
    :param r_drop: float, the section radius of single droplet, not the real radius of the sphere, [dm]
    :return: float, s_drop is the surface area of droplet, [dm2], v_drop is the volume of droplet, [dm3],
    s0 is the section area, [dm2]
    """
    alpha_pom = float(76.8)
    r_real = r_drop / np.sin(alpha_pom)  # [dm]
    height = r_real * (1 - np.cos(alpha_pom))  # [dm]
    s_drop = np.pi * (4 * r_real * height - height ** 2)  # [dm2]
    v_drop = np.pi * height ** 2 * (r_real - height / 3)  # [dm3]
    s0 = np.pi * r_drop ** 2  # [dm2]
    return s_drop, v_drop, s0  # , h_max, s_max, v_max, s1


def watercondense(f=0.1, t=1800, rh=0.8):
    """
    total amount of water condensed with time, according to the partial pressure difference on ground and underground.
    Select Den Haag as the city for average temperature and relative humidity values.
    Calculating the total amount of condensed water, then derive the number of droplets attached on POM section
    :param f: float, flow rate, roughly estimated, [L/s]
    :param t: float, condensation time, assume in half hour the condensation could reach equilibrium with vapourisation.
    :param rh: float, relative humidity, collected online to calculate the water pressure on the ground, [/]
    :return: v_con: float, the total volume of condensed water due to saturated pressure difference and time, [dm3],
    num: float, the number of condensed water droplets, to compute the total amount of mass loss.
    """
    p1 = 25e-3  # saturated water partial pressure at air temperature, [bar]
    p2 = 12.3e-3  # saturated water partial pressure underground, [bar]
    c1 = p1 * 1e2 * rh / (R * (tem + 11))  # water mol on the earth surface, [mol/L]
    c2 = p2 * 1e2 / (R * tem)  # saturated water mol in the pipe, [mol/L]
    c_con = c1 - c2  # water concentration difference for condensation, [mol/L]
    v_con = c_con * f * t * 18 / 1000  # total volume of condensed water in 1800 s, [L]
    num = int(v_con / V_drop)
    return v_con, num


def acid(steps=3600, realtime=3600):
    """
    This function is calculating the time array and accumulated acid concentration, pH values.
    :param steps: float, the time of HCl gas input, default value is 1 hour, [/]
    :param realtime: float, the time for POM degradation process, default value is 1 hour, [/]
    :return: c_water:float, the acid production, without any limitation, [mol/L]
            pH: float, pH values for acid, to indicate when the reaction could happen, [/]
    """
    global c_water, t, t_step, ph
    t_step = 1
    count = 0
    t = np.ones(realtime)
    for i in range(realtime):
        count += 1
        t[i] = count * t_step
    c_water = t * S_drop * abs(N1) / V_drop  # acid concentration in various droplet, [mol/L]
    if c_water.max() > 0.0001:
        if steps < realtime:
            c_water[steps:realtime] = c_water[steps]
    for i in range(t.shape[0]):
        if c_water[i] > 10:
            c_water[i] = 10
    ph = - np.log10(c_water)
    return t, c_water, ph, steps, realtime


def qssa(k11=5e-8, k=3.3e-10, k10=1e-6, mn=30):
    """
    QEA and QSSA assumptions to calculate the reaction rate and acid concentration with consumption, also the mass loss
    with time consuming, the depth of penetration due to POM depolymerization; in this function, the maximum acid
    concentration is considered, assume the HCl gas will not dissolve into the droplet when the acid concentration reach
     10 mol/L, due to maximum HCl acid concentration is around 12 mol/L.
     Another assumption is that the reaction could happen once the pH low to 3.
    :param mn: float, this is the molecular mass of one acetal group, used to estimate the polymer concentration for the
     random law scission.
            k10 = 1e-6  # [L/(mol*s)],
        # from the paper "Influence of Chain Length on the Rate of Hydrolysis of Polyoxymethylene Ethers", the order is 1e-5
        # the random law contains a self-acceleration process, thus the rate constant reduce an order
            k11 = 5e-8  # [1/s], roughly estimating, could be slightly larger
            k = 3.3e-10  # [1/s], according to another paper,
    :return: mass_loss: float, the mass loss for single droplet, [g]
            volume_loss: float, the amount of POM volume loss for single droplet, [dm3]
            acid_cons: float, the acid concentration along time, [mol/L]
            depth: float, the penetration depth along time, [dm]
    """
    density = 1.42 * 1000  # the density of Delrin NC100, [g/dm3]
    cs = density / mn  # reactant groups concentration, [mol/L]
    mass_loss = np.zeros_like(t)
    rs = np.zeros_like(t)
    acid_cons = np.array(c_water)
    for i in range(steps-1):
        # assume the reaction only happen with pH < 3 situation
        if ph[i] < 4:
            rs[i] = ((k10 * cs + k11) / (k11 + k + k10 * cs) - 1) * k10 * cs * acid_cons[i]
            acid_cons[i + 1] = acid_cons[i] - abs(rs[i] * t_step) + S_drop * abs(N1) / V_drop * t_step
            mass_loss[i+1] = mass_loss[i] + abs(rs[i] * t_step * mn)
        else:
            continue
    volume_loss = mass_loss /density
    pH_C = -np.log10(acid_cons)     # the pH value with reaction consumption
    depth = volume_loss / S_rand    # the penetration depth for single droplet, [dm]
    if t.shape[0] > 300000:
        # assume the HCl gas supplying will maintain until the pH is -1.
        try:
            place = np.where(acid_cons > 10)
            loc = int(place[0][0])
        except IndexError:
            loc = steps
        acid_cons[loc:] = acid_cons[loc-1]
        mass_loss[loc:] = mass_loss[loc-1]
        for i in range(loc-1, t.shape[0]-1):
            rs[i] = ((k10 * cs + k11) / (k11 + k + k10 * cs) - 1) * k10 * cs * acid_cons[i]
            acid_cons[i + 1] = acid_cons[i] - abs(rs[i] * t_step)  # + S_drop * abs(N1) / V_drop * t_step
            mass_loss[i + 1] = mass_loss[i] + abs(rs[i] * t_step * mn)
        volume_loss = mass_loss / density
        pH_C = - np.log10(acid_cons)    # the pH value with reaction consumption
        depth = volume_loss / S_rand    # the penetration depth for single droplet
    rs[t.shape[0]-1] = rs[t.shape[0]-2]
    rs = abs(rs)
    return mass_loss, volume_loss, acid_cons, pH_C, depth, rs


# method 1 from page 98, for homogeneous condition, random decomposition with depolymerization
def method_1(kdep=0.002166667, krand=0.000091667):
    """
    :param kdep: the rate constant of depolymerization, [/]
    :param krand: the rate constant of random chain scission, [L/(mol*s)]
    :return: a, float, the mass loss
    """
    a = np.exp((2 * kdep / krand + 2)
               * (np.log(2 - np.exp(-krand * t[0:3600])) - krand * t[0:3600]))

    return a


# method 2 from page 266, reference 92, diffusion-kinetic control
def method_2(kdep=0.002, krand=8e-6, Mn0=70000):
    """
    model 2 from previous research
    :param kdep: rate constant of depolymerization, s^-1
    :param krand: rate constant of random chain scission, s^-1
    :param Mn0:  the
    :return:
    """
    a = np.exp((kdep * t * (1 / Mn0 + krand * c_water / 1e-5)) / (-2.3))
    return a


tem, pre, pwater, phcl = partial_pressure()  # input HCl gas fraction and various operating conditions
x1d, x1_bar, x2d, x2_bar, x3d, x3_bar = average_concentration()
N1 = flux()  # [mol/(dm2*sec)], [dm2/s], [mol/L]
S_drop, V_drop, S_rand = droplet()  # input different droplet radius, in dm^2 and dm^3
V_tot, Number = watercondense()


t0, acid_con0, pH0, steps, realtime0 = acid(900, 2592000)  # input time step and steps, [s], [/]
ans0_mass0, ans0_volume0, ans0_acid0, pH_C0, delta0, rs0= qssa()
ml0 = ans0_mass0 * Number

ans1 = method_1()
ans2 = method_2()

t1, acid_con1, pH1, steps, realtime1 = acid(1800, 2592000)  # input time step and steps, [s], [/]
ans0_mass1, ans0_volume1, ans0_acid1, pH_C1, delta1, rs1 = qssa()
ml1 = ans0_mass1 * Number
t2, acid_con2, pH2, steps, realtime2 = acid(3600, 2592000)  # input time step and steps, [s], [/]
ans0_mass2, ans0_volume2, ans0_acid2, pH_C2, delta2, rs2 = qssa()
ml2 = ans0_mass2 * Number
t3, acid_con3, pH3, steps, realtime3 = acid(14400, 2592000)  # input time step and steps, [s], [/]
ans0_mass3, ans0_volume3, ans0_acid3, pH_C3, delta3, rs3 = qssa()
ml3 = ans0_mass3 * Number
t4, acid_con4, pH4, steps, realtime4 = acid(86400, 2592000)  # input time step and steps, [s], [/]
ans0_mass4, ans0_volume4, ans0_acid4, pH_C4, delta4, rs4 = qssa()
ml4 = ans0_mass4 * Number
t5, acid_con5, pH5, steps, realtime5 = acid(172800, 2592000)  # input time step and steps, [s], [/]
# ans0_mass5, ans0_volume5, ans0_acid5, pH_C5, delta5, rs5 = qssa(5e-8, 1e-8)
# ans0_mass6, ans0_volume6, ans0_acid6, pH_C6, delta6, rs6 = qssa(5e-8, 1e-9)
ans0_mass7, ans0_volume7, ans0_acid7, pH_C7, delta7, rs7 = qssa(5e-8, 3.3e-10)
ml7 = ans0_mass7 * Number

# t6, acid_con6, pH6, steps, realtime6 = acid(172800, 2592000)
# ans0_mass8, ans0_volume8, ans0_acid8, pH_C8, delta8, rs8 = qssa()
# ml8 = ans0_mass8 * Number

'plot part codes'
#
# plt.figure(figsize=(8, 6))
# plt.plot(t0, acid_con0, label="15 min", color="red", linewidth=2)
# plt.plot(t1, acid_con1, label="30 min", color="blue", linewidth=2)
# plt.plot(t2, acid_con2, label="60 min", color="black", linewidth=2)
# plt.plot(t3, acid_con3, label="120 min", color="green", linewidth=2)
# plt.xlabel("Time (min)")
# plt.ylabel("Acid Concentration (mol/L)")
# # plt.title("HCl Acid Concentration with Various HCl Gas Feeding Time")
# plt.xlim(0, realtime0*1.35)
# # plt.ylim(0, acid_con.max()*1.1)
# plt.xticks([0, realtime0/4, realtime0/2, realtime0/4*3, realtime0],
#        ['0', '30', '60', '90', '120'])
# plt.legend(bbox_to_anchor=(0.9, 0.9),
#            bbox_transform=plt.gcf().transFigure)
#
# plt.figure(figsize=(8, 6))
# plt.plot(t3, pH0, label="15 min", color="red", linewidth=6)
# plt.plot(t3, pH_C0, label="15 min with reaction", color="green", linewidth=2)
# plt.plot(t3, pH1, label="30 min", color="blue", linewidth=6)
# plt.plot(t3, pH_C1, label="30 min with reaction", color="yellow", linewidth=2)
# plt.plot(t3, pH2, label="60 min", color="black", linewidth=6)
# plt.plot(t3, pH_C2, label="60 min with reaction", color="white", linewidth=2)
# plt.plot(t3, pH3, label="120 min", color="brown", linewidth=6)
# plt.plot(t3, pH_C3, label="120 min with reaction", color="orange", linewidth=2)
# plt.xlabel("Time (min)")
# plt.ylabel("pH Value")
# # plt.title("HCl Acid Concentration with Various HCl Gas Feeding Time")
# plt.xlim(0, realtime0*1.1)
# # plt.ylim(0, acid_con.max()*1.1)
# plt.xticks([0, realtime0/4, realtime0/2, realtime0/4*3, realtime0],
#        ['0', '30', '60', '90', '120'])
# plt.legend(bbox_to_anchor=(0.48, 0.9),
#            bbox_transform=plt.gcf().transFigure)

# plt.plot(t5, pH5, label="pH level", color="red", linewidth=2)
# # plt.plot(t4, acid_con4, label="7 days", color="blue", linewidth=2)
# # plt.plot(t5, ans0_acid5, label="14 days with degradation", color="green", linewidth=2)
# # plt.plot(t2, acid_con2, label="60 min", color="black", linewidth=2)
# # plt.plot(t3, acid_con3, label="120 min", color="green", linewidth=2)
# plt.xlabel("Time (day)")
# plt.ylabel("pH values")
# # plt.title("HCl Acid Concentration with Various HCl Gas Feeding Time")
# plt.xlim(0, realtime5*1.1)
# plt.ylim(-2, 5)
# plt.xticks([0, realtime5/2, realtime5],
#        ['0', '7', '14'])
# # plt.yticks([-1, 0, +1],
#        # [r'$-1$', r'$0$', r'$+1$'])
# # plt.xticks(np.linspace(0, steps, 6, endpoint=True))
# # plt.yticks(np.linspace(0, acid_con.max()*1.1, 7, endpoint=True))
# plt.legend(bbox_to_anchor=(0.9, 0.9),
#            bbox_transform=plt.gcf().transFigure)


# plt.figure(figsize=(8, 6))
# plt.plot(t0, pH0, label="15 min", color="red", linewidth=2)
# plt.plot(t1, pH1, label="30 min", color="blue", linewidth=2)
# plt.plot(t2, pH2, label="60 min", color="black", linewidth=2)
# plt.plot(t3, pH3, label="120 min", color="green", linewidth=2)
# plt.xlabel("Time (min)")
# plt.ylabel("pH Value")
# # plt.title("pH Accumulation In Water Droplet")
# plt.xlim(0, realtime0*1.1)
# # plt.ylim(-2, 5)
# plt.xticks([0, realtime0/4, realtime0/2, realtime0/4*3, realtime0],
#        ['0', '30', '60', '90', '120'])
# # plt.yticks(np.linspace(-2, 5, 8, endpoint=True))
# plt.legend(bbox_to_anchor=(0.9, 0.9),
#            bbox_transform=plt.gcf().transFigure)

# plt.figure(figsize=(8, 6))
# plt.plot(t5, acid_con5, label="only production", color="red", linewidth=2)
# plt.plot(t5, ans0_acid8, label="hydrolysis process", color="blue", linewidth=2)
# # plt.plot(t5, ans0_acid7, label="1e-9", color="black", linewidth=2)
# # plt.plot(t5, ans0_acid8, label="3.3e-10", color="green", linewidth=2)
# plt.xlabel("Time (day)")
# plt.ylabel("Acid Concentration (mol/L)")
# # plt.title("HCl Acid Concentration with Various HCl Gas Feeding Time")
# plt.xlim(0, realtime5*1.1)
# # plt.ylim(0, acid_con.max()*1.1)
# plt.xticks([0, realtime5/5, realtime5*2/5, realtime5*3/5, realtime5*4/5, realtime5],
#        ['0', '6', '12', '18', '24', '30'])
# # plt.yticks([-1, 0, +1],
# #        [r'$-1$', r'$0$', r'$+1$'])
# plt.legend(bbox_to_anchor=(0.9, 0.9),
#            bbox_transform=plt.gcf().transFigure)


# plt.figure(figsize=(8, 6))
# plt.plot(t5, ans0_acid5, label="1e-7", color="red", linewidth=2)
# plt.plot(t5, ans0_acid6, label="5e-8", color="blue", linewidth=2)
# plt.plot(t5, ans0_acid7, label="1e-9", color="black", linewidth=2)
# # plt.plot(t5, ans0_acid8, label="3.3e-10", color="green", linewidth=2)
# plt.xlabel("Time (day)")
# plt.ylabel("Acid Concentration (mol/L)")
# # plt.title("HCl Acid Concentration with Various HCl Gas Feeding Time")
# plt.xlim(0, realtime5*1.35)
# # plt.ylim(0, acid_con.max()*1.1)
# plt.xticks([0, realtime5/5, realtime5*2/5, realtime5*3/5, realtime5*4/5, realtime5],
#        ['0', '6', '12', '18', '24', '30'])
# # plt.yticks([-1, 0, +1],
# #        [r'$-1$', r'$0$', r'$+1$'])
# plt.legend(bbox_to_anchor=(0.9, 0.85),
#            bbox_transform=plt.gcf().transFigure)

# plt.figure(figsize=(8, 6))
# plt.plot(t5, rs5, label="1e-6", color="red", linewidth=2)
# plt.plot(t5, rs6, label="1e-7", color="blue", linewidth=2)
# plt.plot(t5, rs7, label="5e-8", color="black", linewidth=2)
# plt.plot(t5, rs8, label="1e-9", color="green", linewidth=2)
# plt.xlabel("Time (day)")
# plt.ylabel("Reaction rate (mol/(L*s))")
# # plt.title("HCl Acid Concentration with Various HCl Gas Feeding Time")
# plt.xlim(0, realtime5*1.1)
# # plt.ylim(-1.5, 1)
# plt.xticks([0, realtime5/5, realtime5*2/5, realtime5*3/5, realtime5*4/5, realtime5],
#        ['0', '6', '12', '18', '24', '30'])
# # plt.yticks([-1, 0, +1],
# #        [r'$-1$', r'$0$', r'$+1$'])
# plt.legend(bbox_to_anchor=(0.9, 0.85),
#            bbox_transform=plt.gcf().transFigure)

# plt.figure(figsize=(8, 6))
# plt.plot(t5, rs5, label="1e-8", color="red", linewidth=2)
# plt.plot(t5, rs6, label="1e-9", color="blue", linewidth=2)
# plt.plot(t5, rs7, label="3.3e-10", color="black", linewidth=2)
# plt.plot(t5, rs8, label="1e-11", color="green", linewidth=2)
# plt.xlabel("Time (day)")
# plt.ylabel("Reaction rate (mol/(L*s))")
# # plt.title("HCl Acid Concentration with Various HCl Gas Feeding Time")
# plt.xlim(0, realtime5*1.1)
# # plt.ylim(-1.5, 1)
# plt.xticks([0, realtime5/5, realtime5*2/5, realtime5*3/5, realtime5*4/5, realtime5],
#        ['0', '6', '12', '18', '24', '30'])
# # plt.yticks([-1, 0, +1],
# #        [r'$-1$', r'$0$', r'$+1$'])
# plt.legend(bbox_to_anchor=(0.9, 0.8),
#            bbox_transform=plt.gcf().transFigure)


# plt.figure(figsize=(8, 6))
# plt.plot(t5, pH_C5, label="1e-7", color="red", linewidth=2)
# plt.plot(t5, pH_C6, label="1e-8", color="blue", linewidth=2)
# plt.plot(t5, pH_C7, label="1e-9", color="black", linewidth=2)
# plt.plot(t5, pH_C8, label="3.3e-10", color="green", linewidth=2)
# plt.xlabel("Time (day)")
# plt.ylabel("pH value")
# # plt.title("HCl Acid Concentration with Various HCl Gas Feeding Time")
# plt.xlim(0, realtime5*1.1)
# plt.ylim(-1.5, 1)
# plt.xticks([0, realtime5/5, realtime5*2/5, realtime5*3/5, realtime5*4/5, realtime5],
#        ['0', '6', '12', '18', '24', '30'])
# # plt.yticks([-1, 0, +1],
# #        [r'$-1$', r'$0$', r'$+1$'])
# plt.legend(bbox_to_anchor=(0.9, 0.85),
#            bbox_transform=plt.gcf().transFigure)

# plt.figure(figsize=(8, 6))
# plt.plot(t6, pH6, label="only production", color="red", linewidth=2)
# plt.plot(t6, pH_C8, label="with degradation", color="green", linewidth=2)
# plt.xlabel("Time (day)")
# plt.ylabel("pH value")
# # plt.title("HCl Acid Concentration with Various HCl Gas Feeding Time")
# plt.xlim(0, realtime6*1.1)
# plt.ylim(-1.5, 1)
# plt.xticks([0, realtime6/5, realtime6*2/5, realtime6*3/5, realtime6*4/5, realtime6],
#        ['0', '6', '12', '18', '24', '30'])
# # plt.yticks([-1, 0, +1],
# #        [r'$-1$', r'$0$', r'$+1$'])
# plt.legend(bbox_to_anchor=(0.9, 0.85),
#            bbox_transform=plt.gcf().transFigure)

# #plt.figure(figsize=(8, 6))
# plt.plot(t6, ml8, label="1 month", color="red", linewidth=2)
# # plt.plot(t5, ml1, label="30 min", color="blue", linewidth=2)
# # plt.plot(t5, ml2, label="60 min", color="black", linewidth=2)
# # plt.plot(t5, ml3, label="4 hours", color="brown", linewidth=2)
# # plt.plot(t5, ml4, label="24 hours", color="orange", linewidth=2)
# # plt.plot(t5, ml7, label="48 hours", color="purple", linewidth=2)
# plt.xlabel("Time (day)")
# plt.ylabel("Total Mass loss (g)")
# # plt.title("Mass loss in single droplet")
# plt.xlim(0, realtime6*1.1)
# # plt.ylim(-1.5, 1)
# plt.xticks([0, realtime6/5, realtime6*2/5, realtime6*3/5, realtime6*4/5, realtime6],
#        ['0', '6', '12', '18', '24', '30'])
# # plt.yticks([-1, 0, +1],
# #        [r'$-1$', r'$0$', r'$+1$'])
# plt.legend(bbox_to_anchor=(0.35, 0.9),
#            bbox_transform=plt.gcf().transFigure)

plt.figure(figsize=(8, 6))
plt.plot(t5, delta0, label="15 min", color="red", linewidth=2)
plt.plot(t5, delta1, label="30 min", color="blue", linewidth=2)
plt.plot(t5, delta2, label="60 min", color="black", linewidth=2)
plt.plot(t5, delta3, label="4 hours", color="brown", linewidth=2)
# plt.plot(t5, delta4, label="24 hours", color="orange", linewidth=2)
# plt.plot(t5, delta7, label="48 hours", color="purple", linewidth=2)
plt.xlabel("Time (day)")
plt.ylabel("Degradation Penetration Depth (dm)")
# plt.title("Mass loss in single droplet")
plt.xlim(0, realtime5*1.1)
# plt.ylim(-1.5, 1)
plt.xticks([0, realtime5/5, realtime5*2/5, realtime5*3/5, realtime5*4/5, realtime5],
       ['0', '6', '12', '18', '24', '30'])
# plt.yticks([-1, 0, +1],
#        [r'$-1$', r'$0$', r'$+1$'])
plt.legend(bbox_to_anchor=(0.35, 0.9),
           bbox_transform=plt.gcf().transFigure)
#
# plt.figure(figsize=(8, 6))
# plt.plot(t0, ans1, label="Model 1", color="red", linewidth=2)
# plt.plot(t0, ans2, label="Model 2", color="blue", linewidth=2)
# plt.xlabel("Time(min)")
# plt.ylabel("m/m0")
# plt.xlim(0, realtime0*1.1)
# plt.xticks([0, realtime0/4, realtime0/2, realtime0/4*3, realtime0],
#        ['0', '15', '30', '45', '60'])
# plt.legend(bbox_to_anchor=(0.9, 0.9),
#            bbox_transform=plt.gcf().transFigure)
# plt.draw()


plt.show()

# print h_max, V_drop, V_max, S_drop, S_max, S_rand, delta, ans0_mass
