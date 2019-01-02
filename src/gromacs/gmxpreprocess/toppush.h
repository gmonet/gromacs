/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2018,2019, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

#ifndef GMX_GMXPREPROCESS_TOPPUSH_H
#define GMX_GMXPREPROCESS_TOPPUSH_H

#include "gromacs/fileio/warninp.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/gpp_bond_atomtype.h"

enum class Directive : int;

namespace gmx
{
struct ExclusionBlocks;
} // namespace gmx

void generate_nbparams(int comb, int funct, t_params plist[],
                       gpp_atomtype_t atype,
                       warninp_t wi);

void push_at (struct t_symtab *symtab, gpp_atomtype_t at,
              t_bond_atomtype bat, char *line, int nb_funct,
              t_nbparam ***nbparam, t_nbparam ***pair,
              warninp_t wi);

void push_bt(Directive d, t_params bt[], int nral,
             gpp_atomtype_t at, t_bond_atomtype bat, char *line,
             warninp_t wi);

void push_dihedraltype(Directive d, t_params bt[],
                       t_bond_atomtype bat, char *line,
                       warninp_t wi);

void push_cmaptype(Directive d, t_params bt[], int nral, gpp_atomtype_t at,
                   t_bond_atomtype bat, char *line,
                   warninp_t wi);

void push_nbt(Directive d, t_nbparam **nbt, gpp_atomtype_t atype,
              char *plines, int nb_funct,
              warninp_t wi);

void push_atom(struct t_symtab *symtab,
               t_block         *cgs,
               t_atoms         *at,
               gpp_atomtype_t   atype,
               char            *line,
               int             *lastcg,
               warninp_t        wi);

void push_bond(Directive d, t_params bondtype[], t_params bond[],
               t_atoms *at, gpp_atomtype_t atype, char *line,
               bool bBonded, bool bGenPairs, real fudgeQQ,
               bool bZero, bool *bWarn_copy_A_B,
               warninp_t wi);

void push_cmap(Directive d, t_params bondtype[], t_params bond[],
               t_atoms *at, gpp_atomtype_t atype, char *line,
               warninp_t wi);

void push_vsitesn(Directive d, t_params bond[],
                  t_atoms *at, char *line,
                  warninp_t wi);

void push_mol(int nrmols, t_molinfo mols[], char *pline,
              int *whichmol, int *nrcopies,
              warninp_t wi);

void push_molt(struct t_symtab *symtab, int *nmol, t_molinfo **mol, char *line,
               warninp_t wi);

void push_excl(char *line, gmx::ExclusionBlocks *b2, warninp_t wi);

int add_atomtype_decoupled(struct t_symtab *symtab, gpp_atomtype_t at,
                           t_nbparam ***nbparam, t_nbparam ***pair);
/* Add an atom type with all parameters set to zero (no interactions).
 * Returns the atom type number.
 */

void convert_moltype_couple(t_molinfo *mol, int atomtype_decouple,
                            real fudgeQQ,
                            int couple_lam0, int couple_lam1,
                            bool bCoupleIntra,
                            int nb_funct, t_params *nbp,
                            warninp_t wi);
/* Setup mol such that the B-state has no interaction with the rest
 * of the system, but full interaction with itself.
 */

#endif
