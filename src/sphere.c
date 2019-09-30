/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   sphere.c                                           :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahorker <marvin@42.fr>                     +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2019/09/30 03:59:58 by ahorker           #+#    #+#             */
/*   Updated: 2019/09/30 03:59:58 by ahorker          ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../includes/rt.h"

t_float4	intersect_ray_sphere(t_float4 o, t_float4 d, t_sphere sphere)
{
	float kd[4];
	t_float4 t;
	t_float4 oc;

	t.v[0] = -1;
	t.v[1] = -1;
	t.v[2] = -1;
	t.v[3] = -1;
	oc = gm_vec_subtract(o, sphere.center);
	kd[0] = gm_vec_dot(d, d);
	kd[1] = 2 * gm_vec_dot(oc, d);
	kd[2] = gm_vec_dot(oc, oc) - sphere.data_1.v[1] * sphere.data_1.v[1];
	kd[3] = kd[1] * kd[1] - 4 * kd[0] * kd[2];
	if (kd[3] < 0)
		return (t);
	kd[3] = sqrt(kd[3]);
	t.v[0] = (-kd[1] + kd[3]) / (2 * kd[0]);
	t.v[1] = (-kd[1] - kd[3]) / (2 * kd[0]);
	return (t);
}

t_float4	closest_intersect_sphere(t_rt *rt, t_float4 o,
		t_float4 d, t_float4 lim)
{
	int			iid[2];
	float		closest;
	t_float4	intsect;

	iid[0] = -1;
	iid[1] = -1;
	closest = INF;
	while (++(iid[0]) < rt->obj.count_obj[2])
	{
		if ((int)lim.v[3] == iid[0])
			continue ;
		intsect = intersect_ray_sphere(o, d, rt->obj.spheres[iid[0]]);
		if (intsect.v[0] > lim.v[0] && intsect.v[0] < lim.v[1])
			if (closest > intsect.v[0])
			{
				closest = intsect.v[0];
				iid[1] = iid[0];
			}
		if (intsect.v[1] > lim.v[0] && intsect.v[1] < lim.v[1])
			if (closest > intsect.v[1])
			{
				closest = intsect.v[1];
				iid[1] = iid[0];
			}
	}
	return (gm_init_float(closest, iid[1], 0.f, 0.f));
}

/* ************************************************************************** */
/*                                                                            */
/*           t_mat4 t contains 4 vectors: o - point start ray,                */
/*                    d ray vector, closest and lim.                          */
/*                                                                            */
/*         t_float4 closest contains distance to closest sphere               */
/*                       and id closest sphere.                               */
/*                                                                            */
/*               t_float4 lim contains information about                      */
/*              minimal, maximum intersection calculates                      */
/*             and recursion depth and sphere_intersection                    */
/*          or minus one if this is first ray for this sphere.                */
/*                                                                            */
/* ************************************************************************** */

float	lighting_sphere(t_rt *rt, t_mat4 t, t_sphere s)
{
	float		i;
	float		nl;
	float		rv;
	float		max;
	int 		j;
	t_float4	l;
	t_float4	r;
	t_float4	shadow;

	i = 0.0f;
	j = -1;
	while (++j < rt->obj.count_obj[4])
	{
		if (rt->obj.lights[j].type == ambient)
		{
			i += rt->obj.lights[j].intn;
			continue ;
		}
		else if (rt->obj.lights[j].type == point)
		{
			l = gm_vec_subtract(rt->obj.lights[j].pos, t.v[0]);
			max = 0.f;
		}
		else if (rt->obj.lights[j].type == directional)
		{
			l = rt->obj.lights[j].pos;
			max = INF;
		}
		if (closest_intersect_sphere(rt, t.v[0], l,
				gm_init_float(0.f, max, 0.f, (float)s.data.v[0])).v[1] != -1)
			continue ;
		if (closest_intersect_cone(rt, t.v[0], l,
				gm_init_float(0.f, max, 0.f, -1)).v[1] != -1)
			continue ;
		if (closest_intersect_cylinder(rt, t.v[0], l,
				gm_init_float(0.f, max, 0.f, -1)).v[1] != -1)
			continue ;
		if (closest_intersect_plane(rt, t.v[0], l,
				gm_init_float(0.f, max, 0.f, -1)).v[1] != -1)
			continue ;
		if ((nl = gm_vec_dot(t.v[2], l)) > 0)
			i += rt->obj.lights[j].intn * nl / (gm_vec_len(t.v[2])
					* gm_vec_len(l));
		if (s.data.v[3] > 0)
		{
			r = reflect_ray(l, t.v[2]);
			if ((rv = gm_vec_dot(r, t.v[1])) > 0)
				i += rt->obj.lights[j].intn * pow(rv / (gm_vec_len(r)
						* gm_vec_len(t.v[1])), s.data.v[3]);
		}
	}
	return (i);
}

/* ************************************************************************** */
/*                                                                            */
/*           t_mat4 t contains 4 vectors: o - point start ray,                */
/*                    d ray vector, closest and lim.                          */
/*                                                                            */
/*         t_float4 closest contains distance to closest sphere               */
/*                       and id closest sphere.                               */
/*                                                                            */
/*               t_float4 lim contains information about                      */
/*              minimal, maximum intersection calculates                      */
/*             and recursion depth and sphere_intersection                    */
/*          or minus one if this is first ray for this sphere.                */
/*                                                                            */
/* ************************************************************************** */

t_int4		get_pixel_color_sphere(t_rt *rt, t_mat4 t)
{
	float		light, r;
	int 		color;
	t_int4		argb, local_col;
	t_sphere	s;
	t_float4	n, p;

	s = rt->obj.spheres[(int)t.v[2].v[1]];
	p = gm_vec_add(t.v[0], gm_vec_alp_mult(t.v[1], t.v[2].v[0]));
	n = gm_vec_normalize(gm_vec_subtract(p, s.center));
	light = lighting_sphere(rt, gm_mat_create(p, gm_vec_minus(t.v[1]), n,
			gm_init_float(0, 0, 0, 0)), s);
	color = s.data.v[2];

	argb.v[0] = 0;
	argb.v[1] = (color >> 16 & 0xFF) * light;
	argb.v[2] = (color >> 8 & 0xFF) * light;
	argb.v[3] = (color & 0xFF) * light;

	r = s.data_1.v[0];
	if ((t.v[3].v[2] > 0) && (r > 0))
	{
		t.v[3].v[2] -= 1;
		local_col = trace_ray(rt, gm_mat_create(p, reflect_ray(
				gm_vec_minus(t.v[1]), gm_vec_normalize(
						gm_vec_subtract(p, s.center))),
				gm_init_float(0, 0, 0, 0), t.v[3]));
		argb.v[1] = argb.v[1] * (1 - r) + (local_col.v[1] & 0xFF) * r * 1;
		argb.v[2] = argb.v[2] * (1 - r) + (local_col.v[2] & 0xFF) * r * 1;
		argb.v[3] = argb.v[3] * (1 - r) + (local_col.v[3] & 0xFF) * r * 1;
	}
	return (argb);
}

