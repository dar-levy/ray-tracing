from helper_classes import *
import matplotlib.pyplot as plt

def render_scene(camera, ambient, lights, objects, screen_size, max_depth):
    width, height = screen_size
    ratio = float(width) / height
    screen = (-1, 1 / ratio, 1, -1 / ratio)  # left, top, right, bottom

    image = np.zeros((height, width, 3))

    for i, y in enumerate(np.linspace(screen[1], screen[3], height)):
        for j, x in enumerate(np.linspace(screen[0], screen[2], width)):
            # screen is on origin
            pixel = np.array([x, y, 0])
            origin = camera
            direction = normalize(pixel - origin)
            ray = Ray(origin, direction)

            color = np.zeros(3)
            color = calculate_color(camera, ambient, lights, objects, ray, max_depth, 1)

            # We clip the values between 0 and 1 so all pixel values will make sense.
            image[i, j] = np.clip(color,0,1)

    return image


def calculate_reflected_color(p, camera, ambient, lights, objects, ray, max_depth, level, normal_of_intersection, nearest_object):
    kr = nearest_object.reflection
    reflected_ray = Ray(p, normalize(reflected(ray.direction, normal_of_intersection)))
    return kr * calculate_color(camera, ambient, lights, objects, reflected_ray, max_depth, level)

def calculate_refracted_color(p, camera, ambient, lights, objects, ray, max_depth, level, normal_of_intersection, nearest_object):
    kt = nearest_object.refraction
    if kt > 0:
        n1 = 1.0  # Refractive index of air
        n2 = nearest_object.refractive_index  # Refractive index of the object
        cos_i = np.dot(ray.direction, -normal_of_intersection)
        sin_t2 = (n1 / n2) ** 2 * (1 - cos_i ** 2)
        if sin_t2 > 1:
            # Total internal reflection
            return np.zeros(3)
        cos_t = np.sqrt(1 - sin_t2)
        refracted_direction = (n1 / n2) * ray.direction + (n1 / n2 * cos_i - cos_t) * normal_of_intersection
        refracted_ray_origin = p + normal_of_intersection * (nearest_object.radius + 0.001)
        refracted_ray = Ray(refracted_ray_origin, refracted_direction)
        return kt * calculate_color(camera, ambient, lights, objects, refracted_ray, max_depth, level)
    return np.zeros(3)

def calculate_light_contribution(ambient_color, lights, p, normal_of_intersection, camera, objects, nearest_object, min_distance):
    color = np.float64(ambient_color)
    v = normalize(camera - p)
    for light in lights:
        ray_of_light = light.get_light_ray(p)
        nearest_object_to_light, distance_of_nearease_object_to_light = ray_of_light.nearest_intersected_object(objects)
        if not nearest_object_to_light or distance_of_nearease_object_to_light >= min_distance:
            r = normalize(reflected(ray_of_light.direction, normal_of_intersection))
            diffuse = nearest_object.calc_diffuse(light.get_intensity(p), normal_of_intersection, ray_of_light.direction)
            specular = nearest_object.calc_specular(light.get_intensity(p), v, r)
            color += diffuse + specular

    return color


def calculate_ambient_color(ambient, nearest_object):
    return nearest_object.ambient * ambient


def initialize_ray_trace(ray, nearest_object, min_distance):
    intersection_p = ray.origin + (min_distance * ray.direction)
    normal_of_intersection = nearest_object.compute_normal(intersection_p)
    p = intersection_p + normal_of_intersection / (np.e ** 2)
    return p, normal_of_intersection


def calculate_color(camera, ambient, lights, objects, ray, max_depth, level):
    if level > max_depth:
        return np.zeros(3)
    level += 1

    nearest_object, min_distance = ray.nearest_intersected_object(objects)
    if nearest_object is None:
        return np.zeros(3)

    p, normal_of_intersection = initialize_ray_trace(ray, nearest_object, min_distance)

    ambient_color = calculate_ambient_color(ambient, nearest_object)
    color = calculate_light_contribution(ambient_color, lights, p, normal_of_intersection, camera, objects, nearest_object, min_distance)
    color += calculate_reflected_color(p, camera, ambient, lights, objects, ray, max_depth, level, normal_of_intersection, nearest_object)
    color += calculate_refracted_color(p, camera, ambient, lights, objects, ray, max_depth, level, normal_of_intersection, nearest_object)

    return color


# Write your own objects and lights
# TODO
def your_own_scene():
    """
    An implentation of and 8-ball pool table with a cue ball and 15 solid balls.
    """
    # Cue ball (white)
    cue_ball = Sphere(center=[0, 0, 0.01], radius=0.1)
    cue_ball.set_material([0.9, 0.9, 0.9], [0.9, 0.9, 0.9], [0.8, 0.8, 0.8], 50, 0.2)

    # Solid balls (1-15)
    ball_colors = [[1, 0, 0], 
                   [1, 0.6, 0], [1, 1, 0], 
                   [0, 0, 0], [0.6, 0, 1], [0, 1, 0], 
                   [1, 0, 1], [1, 0.6, 0.2], [0, 0.6, 0], [1, 0.5, 0], 
                   [0.6, 0, 0.6], [1, 0.3, 0.7], [0.4, 0.4, 0], [0.6, 0.6, 0.6], [0.3, 0.3, 0.3]]
    ball_positions = [[0, 0.1, -1], 
                      [0.3, 0.1, -1.6], [-0.3, 0.1, -1.6],
                      [0, 0.1, -2.2], [0.6, 0.1, -2.2], [-0.6, 0.1, -2.2],
                      [0.3, 0.1, -2.8], [-0.3, 0.1, -2.8], [0.9, 0.1, -2.8], [-0.9, 0.1, -2.8],
                      [0, 0.1, -3.4], [0.6, 0.1, -3.4], [-0.6, 0.1, -3.4], [1.2, 0.1, -3.4], [-1.2, 0.1, -3.4]]

    balls = []
    for i in range(len(ball_positions)):
        ball = Sphere(center=ball_positions[i], radius=0.3)
        ball.set_material(ambient=[0.1, 0.1, 0.1], diffuse=ball_colors[i], specular=[0.8, 0.8, 0.8], shininess=50, reflection=0.2)
        balls.append(ball)

    table_surface = Plane([0, 1, 0], [0, -0.3, 0])
    table_surface.set_material([0.2, 0.2, 0.2], [0.2, 0.2, 0.2], [1, 1, 1], 1000, 0.5)


    background = Plane([0, 0, 1], [0, 0, -10])
    background.set_material([0.2, 0.6, 0.2], [0.2, 0.6, 0.2], [0, 0, 0], 100, 0.5)  # Green background

    objects = [table_surface, background] + balls + [cue_ball]

    pointlight = PointLight(intensity=np.array([1, 1, 1]), position=np.array([1, 1.5, 1]),
                            kc=0.1, kl=0.1, kq=0.1)
    spotlight = SpotLight(intensity=np.array([1, 1, 1]), position=np.array([0, 2, -4.5]), direction=np.array([0, -1, 0]), 
                        kc=0.1, kl=0.1, kq=0.1)
    lights = [pointlight, spotlight]

    camera = np.array([0, 1, 1])
    
    return camera, lights, objects

def refracted_scene():
    """
    A scene with a glass sphere and a red ball behind it to the right, showing refraction.
    """
    # Glass sphere
    glass_sphere = Sphere(center=[0, 0, 0], radius=0.3)
    glass_sphere.set_material([0.1, 0.1, 0.1], [0.1, 0.1, 0.1], [0.8, 0.8, 0.8], 50, 0.2, 1, refractive_index=1.52)

    # Red ball
    red_ball = Sphere(center=[1, 0, -5], radius=0.5)
    red_ball.set_material([0.1, 0.1, 0.1], [0.8, 0.1, 0.1], [0.8, 0.8, 0.8], 50, 0.2)

    # Background plane
    background = Plane([0, 0, 1], [0, 0, -10])
    background.set_material([0.8, 0.8, 0.8], [0.8, 0.8, 0.8], [0, 0, 0], 100, 0.5)

    objects = [glass_sphere, red_ball, background]

    # Lights
    pointlight = PointLight(intensity=np.array([1, 1, 1]), position=np.array([-2, 2, 2]), kc=0.1, kl=0.1, kq=0.1)
    lights = [pointlight]

    # Camera
    camera = np.array([0, 0, 2])

    return camera, lights, objects
