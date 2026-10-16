"""Трисектриса Лоншама (варіант 4): r = a / cos 3φ, у декартових: x(x² − 3y²) = a(x² + y²).
Три однакові гілки (симетрія 120°): описуємо гілку k=0, φ ∈ (−30°, 30°), інші — її поворот на k·120°."""
import numpy as np

S3 = np.sqrt(3)
ASYM_DIR = np.radians([30, 90, 150])  # напрями трьох асимптот (cos 3φ = 0)


def rot(pts, k):
    c, s = np.cos(k * 2 * np.pi / 3), np.sin(k * 2 * np.pi / 3)
    return np.asarray(pts) @ np.array([[c, s], [-s, c]])


def branches(a, R, n=1500):
    """Три гілки до відстані R від початку координат (далі — за межами екрана)."""
    pm = np.arccos(min(a / R, 1)) / 3
    p = np.linspace(-pm, pm, n)
    r = a / np.cos(3 * p)
    b = np.c_[r * np.cos(p), r * np.sin(p)]
    return [rot(b, k) for k in range(3)]


def asymptotes(a, L):
    """Відрізки довжини 2L: пряма n·p = d, n = (−sin φ0, cos φ0), d = −a / (3 sin 3φ0)."""
    res = []
    for f in ASYM_DIR:
        d = -a / (3 * np.sin(3 * f))
        p0 = d * np.array([-np.sin(f), np.cos(f)])
        t = np.array([np.cos(f), np.sin(f)])
        res.append(np.array([p0 - L * t, p0 + L * t]))
    return res


def at(a, phi, k=0):
    """Точка кривої, одиничні дотична та нормаль (гілка k, кут φ у радіанах)."""
    c, s = np.cos(3 * phi), np.sin(3 * phi)
    r, dr = a / c, 3 * a * s / c**2
    e, n = np.array([np.cos(phi), np.sin(phi)]), np.array([-np.sin(phi), np.cos(phi)])
    P, T = r * e, dr * e + r * n
    T = T / np.linalg.norm(T)
    return rot(P, k), rot(T, k), rot([-T[1], T[0]], k)


def arc_length(a, phi, n=4000):
    """s = ∫₀^φ a·√(1 + 8 sin²3ψ) / cos²3ψ dψ  (від вершини гілки (a, 0))."""
    p = (np.arange(n) + .5) / n * phi  # метод середніх прямокутників
    return phi / n * np.sum(a * np.sqrt(1 + 8 * np.sin(3 * p) ** 2) / np.cos(3 * p) ** 2)


def sector_area(a, phi):
    """S = ½∫₀^φ r² dψ = a²·tg 3φ / 6 (сектор від вершини до променя φ)."""
    return a * a * np.tan(3 * phi) / 6


def curvature_radius(a, phi):
    """k = −8cos⁴3φ / (a(1 + 8 sin²3φ)^{3/2}) — знак не змінюється, тож точок перегину немає."""
    return a * (1 + 8 * np.sin(3 * phi) ** 2) ** 1.5 / (8 * np.cos(3 * phi) ** 4)


if __name__ == "__main__":
    a = 7.0
    for b in branches(a, 300):  # усі точки лежать на x(x²−3y²) = a(x²+y²)
        x, y = b.T
        assert np.allclose(x * (x**2 - 3 * y**2), a * (x**2 + y**2), rtol=1e-9)
    P, T, N = at(a, 0.3)
    assert abs(T @ N) < 1e-12 and abs(np.linalg.norm(T) - 1) < 1e-12
    e = 1e-5  # дотична — числова похідна
    num = (at(a, 0.3 + e)[0] - at(a, 0.3 - e)[0]) / (2 * e)
    assert np.allclose(num / np.linalg.norm(num), T, atol=1e-6)
    assert np.isclose(curvature_radius(a, 0), a / 8)  # у вершині R = a/8
    h = 1e-4  # довжина малої дуги ≈ хорді
    assert np.isclose(arc_length(a, h), np.linalg.norm(at(a, h)[0] - at(a, 0)[0]), rtol=1e-6)
    psi = (np.arange(100000) + .5) / 100000 * 0.2  # S = ½∫r²dψ чисельно
    assert np.isclose(sector_area(a, 0.2), 0.5 * np.sum((a / np.cos(3 * psi)) ** 2) * 0.2 / 100000, rtol=1e-6)
    for f in ASYM_DIR:  # гілки на нескінченності наближаються до своїх асимптот
        n, d = np.array([-np.sin(f), np.cos(f)]), -a / (3 * np.sin(3 * f))
        dist = []
        for k in range(3):  # φ → f з обох боків, усі гілки: хоч одна точка мусить бути на прямій
            for p in (f - 1e-6, f + 1e-6):
                Q = rot(a / np.cos(3 * p) * np.array([np.cos(p), np.sin(p)]), k)
                dist.append(abs(n @ Q - d))
        assert min(dist) < 1e-3, (f, min(dist))
    assert np.allclose(asymptotes(a, 1)[1][:, 0], -a / 3)  # вертикальна асимптота x = −a/3
    print("ok")
