import numpy as np
from nicegui import ui
import curve as C
import transforms as T

M0 = np.eye(3)  # евклідове перетворення, яке діє лише на криву (сітка й осі нерухомі)
M = M0
view = {'half': 120}
fields = []  # (поле, значення за замовчуванням) — для кнопки «Скинути все»
resetting = False
anim = {'on': False, 'da': 1, 'dphi': 1}


def path(pts):
    return 'M' + ' L'.join(f'{x:.3f} {y:.3f}' for x, y in pts)


def render(*_):
    if resetting:
        return
    h, a = view['half'], p_a.value or 1
    phi, k = np.radians(p_phi.value or 0), int(p_branch.value)
    P, Tg, N = C.at(a, phi, k)
    L = h / 2
    seg = lambda c, col, w, extra='': (f'<path d="{path(T.apply(c, M))}" stroke="{col}" stroke-width="{w}" '
                                       f'fill="none" vector-effect="non-scaling-stroke" {extra}/>')
    out = []
    r = np.arange(-h, h + 1, 10)  # сітка й осі
    grid = ' '.join(path(s) for v in r for s in ([[v, -h], [v, h]], [[-h, v], [h, v]]))
    out.append(f'<path d="{grid}" stroke="#ccc" stroke-width="1" vector-effect="non-scaling-stroke"/>')
    out.append(f'<path d="M{-h} 0 L{h} 0 M0 {-h} L0 {h}" stroke="#333" stroke-width="2" vector-effect="non-scaling-stroke"/>')
    out += [seg(b, '#000', 2) for b in C.branches(a, 3 * h)]
    out.append(f'<circle cx="{T.apply([[0, 0]], M)[0][0]}" cy="{T.apply([[0, 0]], M)[0][1]}" r="{h/120}" fill="#000"/>')  # ізольована точка
    if asym.value:
        out += [seg(s, '#e67e22', 1.5, 'stroke-dasharray="6 4"') for s in C.asymptotes(a, 3 * h)]
    if tang.value:
        out.append(seg([P - L * Tg, P + L * Tg], '#1565c0', 2))
    if norm.value:
        out.append(seg([P - L * N, P + L * N], '#2e7d32', 2))
    (px, py), = T.apply([P], M)
    out.append(f'<circle cx="{px}" cy="{py}" r="{h/80}" fill="#d00"/>')  # поточна точка
    out.append(f'<circle cx="{cx.value or 0}" cy="{cy.value or 0}" r="{h/80}" fill="#b0b"/>')  # центр обертання
    xl, yl = [h * .93, -h * .05], [h * .03, h * .93]
    canvas.set_content(
        f'<svg viewBox="{-h} {-h} {2*h} {2*h}" style="width:100%;height:100vh;background:#fff">'
        f'<g transform="scale(1,-1)">{"".join(out)}</g>'
        f'<g font-size="{h/18}" font-family="Arial" fill="#333">'
        f'<text x="{xl[0]}" y="{-xl[1]}">X</text><text x="{yl[0]}" y="{-yl[1]}">Y</text></g></svg>')
    # числові характеристики (інваріантні відносно евклідових перетворень)
    R = C.curvature_radius(a, phi)
    info['len'].text = f'{C.arc_length(a, phi):.3f}'
    info['area'].text = f'{abs(C.sector_area(a, phi)):.3f}'
    info['R'].text = f'{R:.3f}'
    info['pt'].text = f'({px:.2f}; {py:.2f})'


def apply_m(m):
    global M
    M = M @ m
    render()


def reset():
    global M, resetting
    anim['on'] = False
    M, resetting = M0, True
    for f, v in fields:
        f.value = v
    p_branch.value = 0
    resetting = False
    render()


def tick():
    if not anim['on']:
        return
    global resetting
    lo, hi = sorted((amin.value, amax.value))
    a = p_a.value + anim['da'] * 0.5
    if not lo <= a <= hi:
        anim['da'] *= -1
        a = min(max(a, lo), hi)
    f = p_phi.value + anim['dphi'] * 1.0
    if abs(f) > 25:
        anim['dphi'] *= -1
        f = min(max(f, -25), 25)
    resetting = True  # одне перемальовування замість двох
    p_a.value, p_phi.value = a, f
    resetting = False
    render()


def num(label, value, **kw):
    f = ui.number(label, value=value, **kw).props('dense outlined hide-bottom-space').classes('w-full')
    fields.append((f, value))
    f.on_value_change(render)
    return f


def group(title):
    c = ui.card().props('flat bordered').classes('w-full gap-2')
    with c:
        ui.label(title).classes('font-bold')
    return c


def fgrid(cols, *specs):
    with ui.grid(columns=cols).classes('w-full gap-2'):
        return [num(*s) for s in specs]


def toggle_anim(b):
    anim['on'] = not anim['on']
    b.text = 'Стоп анімація' if anim['on'] else 'Старт анімація'


ui.query('.nicegui-content').classes('p-0 gap-0')
with ui.row().classes('w-full no-wrap gap-0'):
    canvas = ui.html('', sanitize=False).classes('flex-grow min-w-0')
    with ui.column().classes('w-[350px] shrink-0 h-screen overflow-auto p-3 gap-2 bg-grey-2'):
        ui.label('Трисектриса Лоншама: r = a / cos 3φ').classes('font-bold')
        with group('Параметри кривої'):
            p_a, = fgrid(1, ('a', 30))
            p_a.props('min=1 step=1')
        with group('Дослідження точки'):
            p_phi, = fgrid(1, ('φ (°, ±29)', 15))
            p_phi.props('min=-29 max=29 step=1')
            p_branch = ui.toggle([0, 1, 2], value=0, on_change=render).props('dense')
            ui.label('Номер гілки (поворот на 120°)').classes('text-xs')
            tang = ui.checkbox('Дотична', value=True, on_change=render)
            norm = ui.checkbox('Нормаль', value=True, on_change=render)
            asym = ui.checkbox('Асимптоти (3 прямі)', value=True, on_change=render)
            info = {}
            with ui.grid(columns=2).classes('w-full gap-x-2'):
                for key, lbl in [('pt', 'Точка'), ('len', 'Довжина дуги від вершини'),
                                 ('area', 'Площа сектора'), ('R', 'Радіус кривини')]:
                    ui.label(lbl).classes('text-sm')
                    info[key] = ui.label('')
            ui.label('Точок перегину немає (k ≠ 0 всюди). Асимптоти: x = −a/3, x ∓ √3·y = 2a/3. '
                     'Площа між кривою й асимптотою нескінченна.').classes('text-xs')
        with group('Анімація (змінюються a та φ)'):
            amin, amax = fgrid(2, ('a min', 10), ('a max', 60))
            btn = ui.button('Старт анімація', on_click=lambda: toggle_anim(btn)).classes('w-full')
        ui.button('Скинути все', on_click=lambda: (reset(), setattr(btn, 'text', 'Старт анімація'))).classes('w-full')
        with group('Евклідові перетворення'):
            ui.label('Зсув').classes('font-bold')
            dx, dy = fgrid(2, ('dx', 0), ('dy', 0))
            ui.button('Застосувати зсув', on_click=lambda: apply_m(T.shift(dx.value, dy.value))).classes('w-full')
            ui.label('Обертання навколо точки').classes('font-bold')
            ang, cx, cy = fgrid(3, ('Кут (°)', 45), ('Центр X', 40), ('Центр Y', 40))
            ui.button('Застосувати обертання', on_click=lambda: apply_m(T.rotate(cx.value, cy.value, ang.value))).classes('w-full')
        half = ui.number('Вид ±', value=120, min=20, step=10).props('dense outlined hide-bottom-space').classes('w-full')
        half.on_value_change(lambda e: (view.update(half=e.value), render()))

ui.timer(0.05, tick)
render()
ui.run(host='127.0.0.1', port=8080, show=False, reload=False, title='Lab2')
