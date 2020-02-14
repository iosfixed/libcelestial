from jplephem.spk import SPK
from pathlib import Path

root = Path(__file__).parent
kernel_path = root / 'de405.bsp'
kernel = SPK.open(str(kernel_path))


def ephem(JD, target):
    if target == 'moon':
        q_earth, v_earth = kernel[3, 399].compute_and_differentiate(JD)  # Earth relative to Earth Barycenter
        q_moon, v_moon = kernel[3, 301].compute_and_differentiate(JD)  # Moon relative to Earth Barycenter

        q = q_moon - q_earth
        v = (v_moon - v_earth) / 86400  # Convert to km/s

        return q, v

    elif target == 'sun':
        q_sun, v_sun = kernel[0, 10].compute_and_differentiate(JD)  # Sun relative to Solar System Barycenter

        q_earth_bc, v_earth_bc = kernel[0, 3].compute_and_differentiate(JD)  # Earth Barycenter relative to Solar System Barycenter
        q_earth, v_earth = kernel[3, 399].compute_and_differentiate(JD)  # Earth relative to Earth Barycenter

        q = q_earth_bc + q_earth - q_sun
        v = v_earth_bc + v_earth - v_sun

        v /= 86400

        return q, v

    elif target == 'earth':
        q_sun, v_sun = kernel[0, 10].compute_and_differentiate(JD)  # Sun relative to Solar System Barycenter

        q_earth_bc, v_earth_bc = kernel[0, 3].compute_and_differentiate(JD)  # Earth Barycenter relative to Solar System Barycenter
        q_earth, v_earth = kernel[3, 399].compute_and_differentiate(JD)  # Earth relative to Earth Barycenter

        q = q_sun - (q_earth_bc + q_earth)
        v = v_sun - (v_earth_bc + v_earth - v_sun)

        v /= 86400

        return q, v

    raise Exception(f'Unknown target {target}')
