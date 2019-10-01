/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   cone.c                                             :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahorker <marvin@42.fr>                     +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2019/09/30 03:59:58 by ahorker           #+#    #+#             */
/*   Updated: 2019/10/01 22:54:26 by dgorold-         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "../includes/rt.h"

void			get_kd_cone(float *kd, t_float4 co, t_float4 d, t_cone cone)
{
	float		dv;
	float		xv;
	float		dco;

	xv = gm_vec_dot(co, cone.norm);
	dv = gm_vec_dot(d, cone.norm);
	dco = gm_vec_dot(d, co);
	kd[0] = gm_vec_dot(d, d) - dv * dv * (1 + cone.data.v[1] * cone.data.v[1]);
	kd[1] = 2 * (dco - dv * xv * (1 + cone.data.v[1] * cone.data.v[1]));
	kd[2] = gm_vec_dot(co, co) - xv * xv * (1 +
			cone.data.v[1] * cone.data.v[1]);
	kd[3] = kd[1] * kd[1] - 4 * kd[0] * kd[2];
}

t_float4		intersect_ray_cone(t_float4 o, t_float4 d, t_cone cone)
{
	float		kd[4];
	t_float4	t;
	t_float4	co;

	t.v[0] = -1;
	t.v[1] = -1;
	t.v[2] = -1;
	t.v[3] = -1;
	co = gm_vec_subtract(o, cone.center);
	get_kd_cone(&kd[0], co, d, cone);
	if (kd[3] < 0)
		return (t);
	kd[3] = sqrt(kd[3]);
	t.v[0] = (-kd[1] + kd[3]) / (2 * kd[0]);
	t.v[1] = (-kd[1] - kd[3]) / (2 * kd[0]);
	return (t);
}

t_float4		closest_intersect_cone(t_rt *rt, t_float4 o,
		t_float4 d, t_float4 lim)
{
	int			iid[2];
	t_float4	intsect;

	iid[0] = -1;
	iid[1] = -1;
	rt->closest = INF;
	while (++(iid[0]) < rt->obj.count_obj[3])
	{
		if ((int)lim.v[3] == iid[0])
			continue ;
		intsect = intersect_ray_cone(o, d, rt->obj.cones[iid[0]]);
		if (intsect.v[0] > lim.v[0] && intsect.v[0] < lim.v[1])
			if (rt->closest > intsect.v[0])
			{
				rt->closest = intsect.v[0];
				iid[1] = iid[0];
			}
		if (intsect.v[1] > lim.v[0] && intsect.v[1] < lim.v[1])
			if (rt->closest > intsect.v[1])
			{
				rt->closest = intsect.v[1];
				iid[1] = iid[0];
			}
	}
	return (gm_init_float(rt->closest, iid[1], 0.f, 0.f));
}

float			part_one_of_light_cone(t_rt *rt, t_mat4 t, t_cone s, float *k)
{
	if (closest_intersect_cone(rt, t.v[0], rt->l,
			gm_init_float(0.f, k[3], 0.f, (float)s.data.v[0])).v[1] != -1)
		return (0.0f);
	if (closest_intersect_cylinder(rt, t.v[0], rt->l,
			gm_init_float(0.f, k[3], 0.f, -1)).v[1] != -1)
		return (0.0f);
	if (closest_intersect_sphere(rt, t.v[0], rt->l,
			gm_init_float(0.f, k[3], 0.f, -1)).v[1] != -1)
		return (0.0f);
	if (closest_intersect_plane(rt, t.v[0], rt->l,
			gm_init_float(0.f, k[3], 0.f, -1)).v[1] != -1)
		return (0.0f);
	if ((k[1] = gm_vec_dot(t.v[2], rt->l)) > 0)
		k[0] += rt->obj.lights[rt->j].intn * k[1] / (gm_vec_len(t.v[2])
				* gm_vec_len(rt->l));
	if (s.id.v[3] > 0)
	{
		rt->r = reflect_ray(rt->l, t.v[2]);
		if ((k[3] = gm_vec_dot(rt->r, t.v[1])) > 0)
			k[0] += rt->obj.lights[rt->j].intn * pow(k[3] / (gm_vec_len(rt->r)
					* gm_vec_len(t.v[1])), s.id.v[3]);
	}
	return (k[0]);
}

float			lighting_cone(t_rt *rt, t_mat4 t, t_cone s)
{
	float		k[4];

	k[0] = 0.0f;
	k[3] = 0.0f;
	rt->j = -1;
	while (++(rt->j) < rt->obj.count_obj[4])
	{
		if (rt->obj.lights[rt->j].type == ambient)
		{
			k[0] += rt->obj.lights[rt->j].intn;
			continue ;
		}
		else if (rt->obj.lights[rt->j].type == point)
			rt->l = gm_vec_subtract(rt->obj.lights[rt->j].pos, t.v[0]);
		else if (rt->obj.lights[rt->j].type == directional)
		{
			rt->l = rt->obj.lights[rt->j].pos;
			k[3] = INF;
		}
		if ((k[0] = part_one_of_light_cone(rt, t, s, &k[0])) == 0.0f)
			continue ;
	}
	return (k[0]);
}

t_int4			result_color_cone(t_rt *rt, t_mat4 t, t_tmp tmp, t_cone s)
{
	t_int4		l_col;
	float		r;

	r = s.data.v[0];
	if ((t.v[3].v[2] > 0) && (r > 0))
	{
		t.v[3].v[2] -= 1;
		l_col = trace_ray(rt, gm_mat_create(tmp.p, reflect_ray(
				gm_vec_minus(t.v[1]), gm_vec_normalize(
						gm_vec_subtract(tmp.p, s.center))),
								gm_init_float(0, 0, 0, 0), t.v[3]));
		tmp.argb.v[1] = tmp.argb.v[1] * (1 - r) + (l_col.v[1] & 0xFF) * r * 1;
		tmp.argb.v[2] = tmp.argb.v[2] * (1 - r) + (l_col.v[2] & 0xFF) * r * 1;
		tmp.argb.v[3] = tmp.argb.v[3] * (1 - r) + (l_col.v[3] & 0xFF) * r * 1;
	}
	return (tmp.argb);
}

t_int4			get_pixel_color_cone(t_rt *rt, t_mat4 t)
{
	t_tmp		tmp;
	t_cone		s;

	s = rt->obj.cones[(int)t.v[2].v[1]];
	tmp.p = gm_vec_add(t.v[0], gm_vec_alp_mult(t.v[1], t.v[2].v[0]));
	tmp.n = gm_vec_normalize(gm_vec_subtract(
			gm_vec_subtract(tmp.p, s.center), s.norm));
	tmp.light = lighting_cone(rt, gm_mat_create(
			tmp.p, gm_vec_minus(t.v[1]), tmp.n,
			gm_init_float(0, 0, 0, 0)), s);
	tmp.color = s.id.v[2];
	tmp.argb.v[0] = 0;
	tmp.argb.v[1] = (tmp.color >> 16 & 0xFF) * tmp.light;
	tmp.argb.v[2] = (tmp.color >> 8 & 0xFF) * tmp.light;
	tmp.argb.v[3] = (tmp.color & 0xFF) * tmp.light;
	tmp.argb = result_color_cone(rt, t, tmp, s);
	return (tmp.argb);
}
