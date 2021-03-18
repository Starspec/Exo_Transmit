/* This file is part of Exo_Transmit.

    Exo_Transmit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Exo_Transmit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Exo_Transmit.  If not, see <http://www.gnu.org/licenses/>.
*/

/* --------------------------------- atmos.h ------------------------ 

Author: Eliza Kempton (kemptone@grinnell.edu)

--------------------------------------------------------------------- */

#ifndef __ATMOS_H__
#define __ATMOS_H__

/* --- Structure defines atmosphere --------------------------------- */
 
struct Atmos {

  double *lambda, *T, *P, *mu, **kappa;
  double *total, *C, *CH4, *CO, *CO2, *C2H2, *C2H4, *C2H6, *H, 
    *HCN, *H2, *H2CO, *H2O, *H2S, *N, *N2, *NH3, *NO, *NO2, *O, *O2, 
    *O3, *OH, *OCS, *S, *SO2; 

};

#endif /* !__ATMOS_H__ */

/* ------- end ---------------------------- atmos.h ----------------- */
