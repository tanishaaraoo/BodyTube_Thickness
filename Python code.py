import numpy as np
import matplotlib.pyplot as plt

# Constants
pi = 3.14159265359

#center of gravity
cg = 95.6

#center of pressure(distance from tip of nosecone)
cp = 123

#root chord length
root = 22

#total length
Tl= 1660
# Section 1: To find drag force

# Air density (kg/m^3)
p = 1.10

# Vehicle velocity (m/s)
V = 180

# Drag coefficient
Cd = 0.485

# Nosecone base diameter (mm)
D = 182 # mm

# Frontal area of rocket (sq. mm)
A = 0.25 * pi * D**2  # mm^2

# Dynamic pressure (N/mm^2)
q = 0.5 * p * V**2 / 1e6  # N/mm^2

# Drag force (N)
fd = q * Cd * A

# Section 2: To find Normal force on Finset

# Horizontal distance from farthest point of the fin to bodytube (mm)
s = 180

# Length of outer edge (tip chord) (mm)
ct = 60

# Length of inner edge (root chord) (mm)
cr = 220

# Interference coefficient (1 for 3-4 fins, 0.5 for 6 fins)
f = 1

# Radius of bodytube (mm)
R = 89

# Number of fins
N = 4

# Effective angle of attack (radians)
a = 10 * (pi / 180)  # degrees to radians

# Midchord sweep angle (degrees)
b = 0

# Midchord length (mm)
l = 180

# Slope of normal force coefficient at a=0, per radian
Cf = (1 + (f * R) / (s + R)) * (4 * N * (s / D)**2 / (1 + np.sqrt(1 + ((2 * l) / (cr + ct))**2)))

# Normal force acting on the fins (N)
Nfins = q * A * a * Cf

# Section 3: Normal force on nosecone

# Slope of normal force coefficient at a=0, per radian (nosecone)
Cn = 2

#safety factor
SF=1.5

sy= 300000000

#nosecone length

L= 400

# Normal force acting on the nosecone (N)
Nnose = q * A * a * Cn

# Section 4: Distributed Loads

# Distances in mm
Xn = 0.466*L  # Distance from tip of nosecone to Cp of nosecone
X1 = cg-Xn # Distance between Cp of nosecone to Cg of body
X2 = Tl-cg-(root/2)  # Distance between Cg of body to Cp of fins
print(f"Nnose{Nnose:.2f}")
print(f"Nfins{Nfins:.2f}")
print(f"xn:{Xn:.2f}")
print(f"x1:{X1:.2f}")
print(f"x2:{X2:.2f}")
# Distributed loads (N/mm)
w2 = (Nfins * (2 * X2 + X1) - Nnose * X1) / (X2**2 + X2 * X1)
w1 = ((Nnose + Nfins) - (w2 * X2)) / X1

# Section 5: Lateral Shear

x = np.linspace(0, X1 + X2, 5)  # Range of x over which both equations are valid
V = np.zeros_like(x)  # Initialize lateral shear (N)

# Lateral shear for (0 <= x <= X1)
V[x >= 0] = Nnose - w1 * x[x >= 0]
V1 = Nnose - w1 * X1  # Lateral shear at X1

# Lateral shear for (x > X1)
V[x > X1] = V1 - w2 * (x[x > X1] - X1)

# Print all values of V
print("Lateral Shear (V) values:")
for i, val in enumerate(V):
    print(f"x[{i}] = {x[i]:.4f}, V[{i}] = {val:.4f} N")

# Section 6: Bending Moment
M = np.zeros_like(x)  # Initialize bending moment (N-mm)

# Bending moment for (0 <= x <= X1)
M[x >= 0] = Nnose * x[x >= 0] - 0.5 * w1 * x[x >= 0]**2

# Bending moment for (x > X1)
M[x > X1] = V1 * x[x > X1] + w2 * (X1 * x[x > X1] + 0.5 * (X1 + X2)**2 - 0.5 * x[x > X1]**2) - (X1 + X2) * (V1 + w2 * X1)

# Print all values of M
print("\nBending Moment (M) values:")
for i, val in enumerate(M):
    print(f"x[{i}] = {x[i]:.4f}, M[{i}] = {val:.4f} N-mm")

# Section 7: Maximum bending moment
Mmax = np.max(M)
print(f"\nMaximum Bending Moment (Mmax): {Mmax:.4f} N-mm")


# Section modulus (mm^3)
Z = (Mmax*1.5)/3
print(f"Section modulus:{Z:.2f}")

# Finding body thickness (mm)
val = Z*32*D/(pi)
print(f"value: {val}")
d = (D**4 - val)**0.25
thickness = (D - d) / 2

# Display results
print(f"Drag Force: {fd:.2f} N")
print(f"Normal Force on Fins: {Nfins:.2f} N")
print(f"Normal Force on Nosecone: {Nnose:.2f} N")
print(f"Distributed Load w1: {w1:.2f} N/mm")
print(f"Distributed Load w2: {w2:.2f} N/mm")
print(f"Maximum Bending Moment: {Mmax:.2f} N-mm")
print(f"Body Thickness: {thickness:.2f} mm")
print(f"inner diameter:{d} mm")


# Plotting Lateral Shear
plt.figure(1)
plt.plot(x, V)
plt.xlabel('x (mm)')
plt.ylabel('Lateral Shear (N)')
plt.title('Lateral Shear vs x')
plt.grid(True)

# Plotting Bending Moment
plt.figure(2)
plt.plot(x, M)
plt.xlabel('x (mm)')
plt.ylabel('Bending Moment (N-mm)')
plt.title('Bending Moment vs x')
plt.grid(True)

plt.show()
