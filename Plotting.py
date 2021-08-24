from numpy import array, zeros
from matplotlib.pyplot import figure, plot, yscale, xlabel, ylabel, title, legend, savefig, show

n = 0
with open("cosmology.txt", "r") as infile:
    for line in infile:
        n += 1

x = zeros(n)
eta = zeros(n)
H = zeros(n)
Hp = zeros(n)
dHp = zeros(n)
ddHp = zeros(n)
OB = zeros(n)
OCDM = zeros(n)
OL = zeros(n)
OR = zeros(n)
ON = zeros(n)
OK = zeros(n)
t = zeros(n)

with open("cosmology.txt", "r") as infile:
    i = 0
    for line in infile:
        x[i], eta[i], H[i], Hp[i], dHp[i], ddHp[i], OB[i], OCDM[i], OL[i], OR[i], ON[i], OK[i], t[i] = line.split()
        i += 1

print(eta * Hp / (3 * 10 ** 8))

figure()
plot(t / (60 * 60 * 24 * 365.25 * 1e9), x)
xlabel("t [Gy]")
ylabel(r"$\ln(a)$")
savefig("Figure01.png")

figure()
plot(x, eta / (3 * 10 ** 22))
xlabel(r"$\ln(a)$")
ylabel(r"Conformal time $\eta$ [Mpc]")
yscale("log")
savefig("Figure02.png")

figure()
plot(x, H)
yscale("log")
xlabel(r"$\ln(a)$")
ylabel(r"Hubble parameter [$s^{-1}$]")
savefig("Figure03.png")

figure()
plot(x, Hp)
yscale("log")
xlabel(r"$\ln(a)$")
ylabel(r"Scaled Hubble parameter $H'$")
savefig("Figure04.png")

figure()
plot(x, - 1 / Hp * dHp)
xlabel(r"$\ln(a)$")
ylabel(r"$- \frac{1}{H'} \frac{dH'}{dx}$")
yscale("log")
savefig("Figure05.png")

figure()
plot(x, 1 / Hp * ddHp)
xlabel(r"$\ln(a)$")
ylabel(r"$\frac{1}{H'} \frac{d^2H'}{dx^2}$")
yscale("log")
savefig("Figure06.png")

figure()
plot(x, OCDM + OB)
plot(x, OL)
plot(x, OR + ON)
legend([r"$\Omega_{CDM} + \Omega_{b}$", r"$\Omega_{\Lambda}$", r"$\Omega_{\gamma} + \Omega_{\nu}$"])
xlabel(r"$\ln(a)$")
ylabel(r"$\Omega_X$")
savefig("Figure07.png")

show()