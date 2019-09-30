/* ************************************************************************** */
/*                                                                            */
/*                                                        :::      ::::::::   */
/*   gm_vec_mult.c                                      :+:      :+:    :+:   */
/*                                                    +:+ +:+         +:+     */
/*   By: ahorker <marvin@42.fr>                     +#+  +:+       +#+        */
/*                                                +#+#+#+#+#+   +#+           */
/*   Created: 2019/05/02 21:12:10 by ahorker           #+#    #+#             */
/*   Updated: 2019/05/02 21:12:10 by ahorker          ###   ########.fr       */
/*                                                                            */
/* ************************************************************************** */

#include "include/libgm.h"

t_float4	gm_vec_mult(t_float4 v1, t_float4 v2)
{
	t_float4	v;

	v.v[0] = v1.v[1] * v2.v[2] - v1.v[2] * v2.v[1] + 0.0;
	v.v[1] = v1.v[2] * v2.v[0] - v1.v[0] * v2.v[2] + 0.0;
	v.v[2] = v1.v[0] * v2.v[1] - v1.v[1] * v2.v[0] + 0.0;
	v.v[3] = 0;
	return (v);
}