/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   rt.h                                               :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahorker <ahorker@student.42.fr>            +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2018/11/20 22:19:38 by ahorker           #+#    #+#             */
/*   Updated: 2019/09/30 22:10:21 by dgorold-         ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#ifndef RT_H
# define RT_H

# include "../libft/includes/libft.h"
# include "../libgm/include/libgm.h"

# include <mlx.h>
# include <math.h>
# include <fcntl.h>
# include <stdio.h>
# include <pthread.h>

//# define SC_X		(1080)
//# define SC_Y		(1080)
//# define SC_X_2		(540)
//# define SC_Y_2		(540)
//# define VW			(1)
//# define VH			(1)
//# define D			(1)
# define INF		(10000000.0f)
# define BG_COLOR	(0xFFFFFF)
# define REC_DEPTH	(20.f)

typedef struct		s_mlx
{
	void			*mlx;
	void			*win;
	void			*img;
	int				bit_per_pix;
	int				size_len;
	int				endian;
	char			*string;
}					t_mlx;

typedef struct		s_triangle
{
	t_int4			id;
	t_float4		point[3];
	t_float4		norm;
	int				color;
}					t_triangle;

typedef struct		s_plane
{
	t_int4			id;
	t_float4		point;
	t_float4		norm;
	t_float4		data;
}					t_plane;


typedef enum		e_light_type
{					ambient, point, directional
}					t_light_type;

typedef enum		e_obj_type
{					cone, sphere, plane, cylinder
}					t_obj_type;


typedef struct		s_light
{
	t_int4			id;
	t_light_type	type;
	t_float4		pos;
	float			intn;
}					t_light;

/* ************************************************************************** */
/*                 t_int4 data contains information about                     */
/*                    id, global_id, color and specular                       */
/*                              for this sphere.                              */
/*               t_float4 data_1 contains information about                   */
/*                           reflection and radius				              */
/*                              for this sphere.                              */
/* ************************************************************************** */

typedef struct		s_sphere
{
	t_int4			data;
	t_float4		data_1;
	t_float4		center;
}					t_sphere;

/* ************************************************************************** */
/*                  t_int4 id contains information about                      */
/*                    id, global_id, color and specular                       */
/*                            for this cone.                                  */
/*                t_float4 data contains information about                    */
/*                 reflection and tan_tetta for this cone.                    */
/* ************************************************************************** */

typedef struct		s_cone
{
	t_int4			id;
	t_float4		center;
	t_float4		norm;
	t_float4		data;
}					t_cone;


/* ************************************************************************** */
/*                  t_int4 id contains information about                      */
/*                    id, global_id, color and specular                       */
/*                            for this cylinder.                              */
/*                   t_float4 data contains information                       */
/*             about reflection and radius for this cylinder.                 */
/* ************************************************************************** */

typedef struct		s_cylinder
{
	t_int4			id;
	t_float4		center;
	t_float4		norm;
	t_float4		data;
}					t_cylinder;

typedef struct		s_obj
{
	int 			count_obj[5];
	t_plane			planes[6];
	t_cylinder		cylinders[3];
	t_sphere		spheres[4];
	t_cone			cones[3];
	t_light			lights[5];
}					t_obj;

typedef struct		s_tmp
{

}					t_tmp;

typedef struct		s_rt
{
	t_mlx			mlx;
	t_int4			**image;
	t_obj			obj;
	t_int4			scr;
	t_float4		cnv;
	t_mat4			t;
	t_float4		o;
}					t_rt;


void	draw_point(t_rt *rt, int x, int y, t_int4 p);
t_float4	rotation(t_float4 main_axis, t_float4 rot_axis, float angle);

void	init_mlx(t_rt *rt);
void	image_to_win(t_rt *rt);
int		deal_key(int key, t_rt *param);

void	print_vect(t_float4 x);
void	print_vect_int(t_int4 x);
t_float4	reflect_ray(t_float4 v, t_float4 n);

t_int4	trace_ray(t_rt *rt, t_mat4 t);
void	calc(t_rt *rt);

t_light		init_light(t_int4 id, t_light_type type, t_float4 vec, float intens);
t_cone		init_cone(t_int4 id, t_float4 norm, t_float4 center, t_float4 data);
t_sphere	init_sphere(t_int4 data, t_float4 data_1, t_float4 center);
t_cylinder	init_cylinder(t_int4 id, t_float4 norm, t_float4 center, t_float4 data);
t_plane		init_plane(t_int4 id, t_float4 point, t_float4 norm, t_float4 data);

t_float4	closest_intersect_cone(t_rt *rt, t_float4 o, t_float4 d, t_float4 lim);
t_float4	closest_intersect_sphere(t_rt *rt, t_float4 o, t_float4 d, t_float4 lim);
t_float4	closest_intersect_cylinder(t_rt *rt, t_float4 o, t_float4 d, t_float4 lim);
t_float4	closest_intersect_plane(t_rt *rt, t_float4 o, t_float4 d, t_float4 lim);

t_int4		get_pixel_color_cone(t_rt *rt, t_mat4 t);
t_int4		get_pixel_color_sphere(t_rt *rt, t_mat4 t);
t_int4		get_pixel_color_cylinder(t_rt *rt, t_mat4 t);
t_int4		get_pixel_color_plane(t_rt *rt, t_mat4 t);


#endif
