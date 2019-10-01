/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   main.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: dgorold- <dgorold-@student.42.fr>          +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2019/09/30 23:17:28 by dgorold-          #+#    #+#             */
/*   Updated: 2019/10/01 20:15:48 by dgorold-         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "rt.h"

void			draw_point(t_rt *rt, int x, int y, t_int4 p)
{
	int			ind;

	if (y > rt->scr.v[1] - 1 || y <= 0 || x > rt->scr.v[0] - 1 || x <= 0)
		return ;
	ind = (rt->mlx.bit_per_pix / 8) * x + rt->mlx.size_len * y;
	rt->mlx.string[ind] = p.v[3] & 0xFF;
	rt->mlx.string[ind + 1] = p.v[2] & 0xFF;
	rt->mlx.string[ind + 2] = p.v[1] & 0xFF;
	rt->mlx.string[ind + 3] = p.v[0] & 0xFF;
}

void			init_mlx(t_rt *rt)
{
	rt->mlx.mlx = mlx_init();
	rt->mlx.win = mlx_new_window(rt->mlx.mlx, rt->scr.v[0], rt->scr.v[1], "RT");
	rt->mlx.img = mlx_new_image(rt->mlx.mlx, rt->scr.v[0], rt->scr.v[1]);
	rt->mlx.string = mlx_get_data_addr(rt->mlx.img,
			&(rt->mlx.bit_per_pix), &(rt->mlx.size_len), &(rt->mlx.endian));
}

t_float4		reflect_ray(t_float4 v, t_float4 n)
{
	return (gm_vec_subtract(gm_vec_alp_mult(n, 2 * gm_vec_dot(n, v)), v));
}

t_obj_type		switch_obj(t_rt *rt, t_obj_type x, t_mat4 *t)
{
	t_float4	t1;

	t->v[2] = closest_intersect_sphere(rt, t->v[0], t->v[1], t->v[3]);
	t1 = closest_intersect_cone(rt, t->v[0], t->v[1], t->v[3]);
	if (t->v[2].v[0] > t1.v[0])
	{
		x = cone;
		t->v[2] = t1;
	}
	t1 = closest_intersect_cylinder(rt, t->v[0], t->v[1], t->v[3]);
	if (t->v[2].v[0] > t1.v[0])
	{
		x = cylinder;
		t->v[2] = t1;
	}
	t1 = closest_intersect_plane(rt, t->v[0], t->v[1], t->v[3]);
	if (t->v[2].v[0] > t1.v[0])
	{
		x = plane;
		t->v[2] = t1;
	}
	return (x);
}

t_int4			trace_ray(t_rt *rt, t_mat4 t)
{
	t_obj_type	x;

	x = sphere;
	x = switch_obj(rt, x, &t);
	if (t.v[2].v[1] == -1)
		return (gm_init_int(BG_COLOR >> 24 & 0xFF, BG_COLOR >> 16 & 0xFF,
				BG_COLOR >> 8 & 0xFF, BG_COLOR & 0xFF));
	else if (x == sphere)
		return (get_pixel_color_sphere(rt, t));
	else if (x == cone)
		return (get_pixel_color_cone(rt, t));
	else if (x == cylinder)
		return (get_pixel_color_cylinder(rt, t));
	else if (x == plane)
		return (get_pixel_color_plane(rt, t));
	return (gm_init_int(0xFF, 0xFF, 0xFF, 0xFF));
}

t_sphere		init_sphere(t_int4 data, t_float4 data_1, t_float4 center)
{
	t_sphere	a;

	a.center = center;
	a.data = data;
	a.data_1 = data_1;
	return (a);
}

t_cone			init_cone(t_int4 id, t_float4 center,
		t_float4 norm, t_float4 data)
{
	t_cone		a;

	a.center = center;
	a.norm = norm;
	a.id = id;
	a.data = data;
	return (a);
}

t_cylinder		init_cylinder(t_int4 id, t_float4 center,
		t_float4 norm, t_float4 data)
{
	t_cylinder	a;

	a.center = center;
	a.norm = norm;
	a.id = id;
	a.data = data;
	return (a);
}

t_plane			init_plane(t_int4 id, t_float4 point,
		t_float4 norm, t_float4 data)
{
	t_plane		a;

	a.point = point;
	a.norm = norm;
	a.id = id;
	a.data = data;
	return (a);
}

t_light			init_light(t_int4 id, t_light_type type,
		t_float4 vec, float intns)
{
	t_light		a;

	a.id = id;
	a.type = type;
	a.pos = vec;
	a.intn = intns;
	return (a);
}

void			calc(t_rt *rt)
{
	int			x;
	int			y;
	t_float4	d;
	t_int4		color;

	x = -rt->scr.v[2] - 1;
	while (++x < rt->scr.v[2] + 1)
	{
		y = -rt->scr.v[3] - 1;
		while (++y < rt->scr.v[3] + 1)
		{
			d = gm_vec_subtract(gm_init_float(
					(float)(x * (rt->cnv.v[0] / rt->scr.v[0]) + rt->o.v[0]),
					(float)(y * (rt->cnv.v[1] / rt->scr.v[1]) + rt->o.v[1]),
					rt->cnv.v[2] + rt->o.v[2], 0), rt->o);
			d = gm_mat_pnt_mult(d, rt->t);
			color = trace_ray(rt, gm_mat_create(rt->o, d,
					gm_init_float(0, 0, 0, 0),
					gm_init_float(1, INF, REC_DEPTH, -1)));
			draw_point(rt, rt->scr.v[2] + x,
					rt->scr.v[3] - y, color);
		}
	}
}

void			image_to_win(t_rt *rt)
{
	rt->mlx.img = mlx_new_image(rt->mlx.mlx, rt->scr.v[0], rt->scr.v[1]);
	rt->mlx.string = mlx_get_data_addr(rt->mlx.img,
			&(rt->mlx.bit_per_pix), &(rt->mlx.size_len), &(rt->mlx.endian));
	calc(rt);
	mlx_clear_window(rt->mlx.mlx, rt->mlx.win);
	mlx_put_image_to_window(rt->mlx.mlx, rt->mlx.win, rt->mlx.img, 0, 0);
	mlx_destroy_image(rt->mlx.mlx, rt->mlx.img);
}

t_float4		rotation(t_float4 ax_1, t_float4 ax_2, float angle)
{
	double		tmp[6];

	tmp[3] = cos(angle);
	tmp[4] = 1 - tmp[3];
	tmp[5] = sin(angle);
	tmp[0] = ax_2.v[0] * (tmp[3] + tmp[4] * ax_1.v[0] * ax_1.v[0])
			+ ax_2.v[1] * (tmp[4] * ax_1.v[0] * ax_1.v[1] - tmp[5] * ax_1.v[2])
			+ ax_2.v[2] * (tmp[4] * ax_1.v[0] * ax_1.v[2] + tmp[5] * ax_1.v[1]);
	tmp[1] = ax_2.v[0] * (tmp[4] * ax_1.v[0] * ax_1.v[1] + tmp[5] * ax_1.v[2])
			+ ax_2.v[1] * (tmp[3] + tmp[4] * ax_1.v[1] * ax_1.v[1])
			+ ax_2.v[2] * (tmp[4] * ax_1.v[1] * ax_1.v[2] - tmp[5] * ax_1.v[0]);
	tmp[2] = ax_2.v[0] * (tmp[4] * ax_1.v[0] * ax_1.v[2] - tmp[5] * ax_1.v[1])
			+ ax_2.v[1] * (tmp[4] * ax_1.v[1] * ax_1.v[2] + tmp[5] * ax_1.v[0])
			+ ax_2.v[2] * (tmp[3] + tmp[4] * ax_1.v[2] * ax_1.v[2]);
	ax_2.v[0] = tmp[0];
	ax_2.v[1] = tmp[1];
	ax_2.v[2] = tmp[2];
	return (ax_2);
}

void			left_right_down_up_deal_key(int key, t_rt *rt)
{
	if (key == 123)
	{
		rt->t.v[0] = rotation(rt->t.v[1], rt->t.v[0], 0.15708);
		rt->t.v[2] = rotation(rt->t.v[1], rt->t.v[2], 0.15708);
	}
	else if (key == 124)
	{
		rt->t.v[0] = rotation(rt->t.v[1], rt->t.v[0], -0.15708);
		rt->t.v[2] = rotation(rt->t.v[1], rt->t.v[2], -0.15708);
	}
	else if (key == 125)
	{
		rt->t.v[2] = rotation(rt->t.v[0], rt->t.v[2], 0.15708);
		rt->t.v[1] = rotation(rt->t.v[0], rt->t.v[1], 0.15708);
	}
	else if (key == 126)
	{
		rt->t.v[2] = rotation(rt->t.v[0], rt->t.v[2], -0.15708);
		rt->t.v[1] = rotation(rt->t.v[0], rt->t.v[1], -0.15708);
	}
}

void			letters_deal_key(int key, t_rt *rt)
{
	if (key == 13)
		rt->o.v[2] += 1;
	else if (key == 1)
		rt->o.v[2] -= 1;
	else if (key == 2)
		rt->o.v[0] += 1;
	else if (key == 0)
		rt->o.v[0] -= 1;
	else if (key == 116)
		rt->o.v[1] += 1;
	else if (key == 121)
		rt->o.v[1] -= 1;
}

int				deal_key(int key, t_rt *rt)
{
	if (key == 53)
		exit(0);
	letters_deal_key(key, rt);
	if (key == 12)
	{
		rt->t.v[0] = rotation(rt->t.v[2], rt->t.v[0], 0.15708);
		rt->t.v[1] = rotation(rt->t.v[2], rt->t.v[1], 0.15708);
	}
	else if (key == 14)
	{
		rt->t.v[0] = rotation(rt->t.v[2], rt->t.v[0], -0.15708);
		rt->t.v[1] = rotation(rt->t.v[2], rt->t.v[1], -0.15708);
	}
	left_right_down_up_deal_key(key, rt);
	image_to_win(rt);
	return (0);
}

void			first_part_of_the_main(t_rt *rt)
{
	rt->t = gm_mat_create(
			gm_init_float(1, 0, 0, 0),
			gm_init_float(0, 1, 0, 0),
			gm_init_float(0, 0, 1, 0),
			gm_init_float(0, 0, 0, 1));
	rt->obj.spheres[0] = init_sphere(gm_init_int(0, 0, 0xFF0000, 1000),
			gm_init_float(0., 1, 0, 0),
			gm_init_float(0, 0, 3, 0));
	rt->obj.count_obj[2]++;
	rt->obj.spheres[1] = init_sphere(gm_init_int(1, 1, 0x0000FF, 500),
			gm_init_float(.3, 1, 0, 0),
			gm_init_float(2, 2, 3, 0));
	rt->obj.count_obj[2]++;
	rt->obj.spheres[2] = init_sphere(gm_init_int(2, 2, 0x00FF00, 10),
			gm_init_float(.0, 1, 0, 0),
			gm_init_float(-2, 2, 3, 0));
	rt->obj.count_obj[2]++;
	rt->obj.planes[0] = init_plane(gm_init_int(0, 3, 0xFFFFFF, 100),
			gm_init_float(0., -0.05, 0, 0),
			gm_init_float(0, 1, 0, 0),
			gm_init_float(0.9, 0, 0, 0));
}

void			second_part_of_the_main(t_rt *rt)
{
	rt->obj.count_obj[0]++;
	rt->obj.lights[4] = init_light(gm_init_int(0, 4, 0, 0), ambient,
			gm_init_float(0, 0, 0, 0), 0.2);
	rt->obj.count_obj[4]++;
	rt->obj.lights[0] = init_light(gm_init_int(1, 5, 0, 0), point,
			gm_init_float(0.3, 1.8, 2.5, 0), .2);
	rt->obj.count_obj[4]++;
	rt->obj.lights[1] = init_light(gm_init_int(2, 6, 0, 0), point,
			gm_init_float(-2, 5, 2, 0), 0.2);
	rt->obj.count_obj[4]++;
	rt->obj.lights[2] = init_light(gm_init_int(3, 7, 0, 0), directional,
			gm_init_float(1, 4, 4, 0), 0.1);
	rt->obj.count_obj[4]++;
	rt->obj.lights[3] = init_light(gm_init_int(4, 8, 0, 0), directional,
			gm_init_float(-1, 4, 4, 0), 0.2);
	rt->obj.count_obj[4]++;
	rt->obj.cones[0] = init_cone(gm_init_int(0, 9, 0xFF00FF, 100),
			gm_init_float(0, 2, 4, 0),
			gm_init_float(0, -1, 0, 0),
			gm_init_float(0.2, 0.26794919243, 0, 0));
	rt->obj.count_obj[3]++;
}

void			third_part_of_the_main(t_rt *rt)
{
	rt->obj.cylinders[0] = init_cylinder(gm_init_int(0, 10, 0x00FFFF, 1000),
			gm_init_float(-3, 2, 4, 0),
			gm_init_float(0., -1, 0, 0),
			gm_init_float(0.3, 1, 0, 0));
	rt->obj.count_obj[1]++;
	rt->obj.planes[1] = init_plane(gm_init_int(1, 11, 0xFFFFFF, 1000),
			gm_init_float(0., 0., 7, 0),
			gm_init_float(0., 0, -1, 0),
			gm_init_float(0.9, 0., 0., 0.));
	rt->obj.count_obj[0]++;
	rt->obj.spheres[3] = init_sphere(gm_init_int(3, 12, 0xFFFFFF, 10),
			gm_init_float(.0, 0.1, 0, 0),
			gm_init_float(0, 0, 0, 0));
	rt->obj.count_obj[2]++;
	rt->obj.planes[2] = init_plane(gm_init_int(2, 13, 0xFFFFFF, 1000),
			gm_init_float(0., 7, 0, 0),
			gm_init_float(0., -1, 0, 0),
			gm_init_float(0.9, 0., 0., 0.));
	rt->obj.count_obj[0]++;
	rt->obj.planes[3] = init_plane(gm_init_int(3, 14, 0x000000, 1000),
			gm_init_float(10, 0., 0, 0),
			gm_init_float(-1, 0, 0, 0),
			gm_init_float(0.9, 0., 0., 0.));
}

void			fourth_part_of_the_main(t_rt *rt)
{
	rt->obj.count_obj[0]++;
	rt->obj.planes[4] = init_plane(gm_init_int(4, 15, 0x000000, 1000),
			gm_init_float(-10., 0., 0, 0),
			gm_init_float(1, 0, 0, 0),
			gm_init_float(0.9, 0., 0., 0.));
	rt->obj.count_obj[0]++;
	rt->obj.planes[5] = init_plane(gm_init_int(4, 15, 0xFFFFFF, 1000),
			gm_init_float(0, 0., -15, 0),
			gm_init_float(0, 0, -1, 0),
			gm_init_float(0.9, 0., 0., 0.));
	rt->obj.count_obj[0]++;
}

int				main(int argc, char *argv[])
{
	t_rt		rt;

	rt.o = gm_init_float(.0f, 3.0f, -10.f, .0f);
	ft_bzero(rt.obj.count_obj, 5 * sizeof(int));
	if (argc == 1)
	{
		first_part_of_the_main(&rt);
		second_part_of_the_main(&rt);
		third_part_of_the_main(&rt);
		fourth_part_of_the_main(&rt);
		rt.cnv = gm_init_float(0.001f * 16.0f / 9.0f, 0.001f, 0.001f, 0);
		rt.scr = gm_init_int(320, 180, 160, 90);
		init_mlx(&rt);
		calc(&rt);
		image_to_win(&rt);
		mlx_hook(rt.mlx.win, 2, 5, deal_key, &rt);
		mlx_loop(rt.mlx.mlx);
	}
	return (0);
}
