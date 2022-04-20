/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017,2018,2019, The GROMACS development team.
 * Copyright (c) 2020,2021, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Defines data structure and utilities for extforce fields
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "extforcefield.h"

#include <cmath>

#include <memory>
#include <string>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdoutputprovider.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/mdtypes/mdatom.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainerwithsections.h"
#include "gromacs/options/optionsection.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/keyvaluetreebuilder.h"
#include "gromacs/utility/keyvaluetreetransform.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/strconvert.h"

namespace gmx
{

namespace
{

/*! \internal
 * \brief Describes an applied extforce field in a coordinate
 * dimension.
 *
 * Can compute the applied field strength at a time, and supports
 * operations to read and form corresponding .mdp contents.
 */
class ExtForceFieldDimension
{
public:
    /*! \brief
     * Adds an option section to specify parameters for this field component.
     */
    void initMdpOptions(IOptionsContainerWithSections* options, const char* sectionName)
    {
        auto section = options->addSection(OptionSection(sectionName));
        section.addOption(RealOption("E0").store(&a_));
        section.addOption(RealOption("omega").store(&omega_));
        section.addOption(RealOption("t0").store(&t0_));
        section.addOption(RealOption("sigma").store(&sigma_));
        section.addOption(RealOption("zmin").store(&zmin_));
        section.addOption(RealOption("zmax").store(&zmax_));
        if (sigma_ <= 0 && t0_ != 0.0)
        {
            GMX_THROW(
                    InvalidInputError("Non-pulsed field (sigma = 0) ignores the value of t0. "
                                      "Please, set t0 to 0 to avoid this error."));
        }
        if (zmin_ > zmax_)
        {
            GMX_THROW(
                    InvalidInputError("zmin must be greater than zmax."));
        }
    }
    /*! \brief
     * Creates mdp parameters for this field component.
     */
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder, const std::string& name) const
    {
        builder->addUniformArray<real>("extforce-field-" + name, { a_, omega_, t0_, sigma_, zmin_, zmax_ });
    }

    /*! \brief Evaluates this field component at given time.
     *
     * \param[in] t The time to evualate at
     * \return The extforce field
     */
    real evaluate(real t) const
    {
        if (sigma_ > 0)
        {
            return a_ * (std::cos(omega_ * (t - t0_)) * std::exp(-square(t - t0_) / (2.0 * square(sigma_))));
        }
        else if (omega_ > 0){
            return a_ * std::cos(omega_ * t);
        }
        else
        {
            return a_;
        }
    }

    
    //! Return the amplitude
    real a() const { return a_; }

    //! Return the lower value of z
    real zmin() const { return zmin_; }

    //! Return the higher value of z
    real zmax() const { return zmax_; }

private:
    //! Coeffient (V / nm)
    real a_ = 0;
    //! Frequency (1/ps)
    real omega_ = 0;
    //! Central time point (ps) for pulse
    real t0_ = 0;
    //! Width of pulse (ps, if zero there is no pulse)
    real sigma_ = 0;
    //! Lower value of z coordinate (nm) for which the field is applied.
    real zmin_ = 0;
    //! Highest value of z coordinate (nm) for which the field is applied.
    real zmax_ = 0;
};

/*! \internal
 * \brief Describe time dependent extforce field
 *
 * Class that implements a force to be evaluated in mdrun.
 * The extforce field can be pulsed and oscillating, simply
 * oscillating, or static, in each of X,Y,Z directions.
 */
class ExtForceField final : public IMDModule, public IMdpOptionProvider, public IMDOutputProvider, public IForceProvider
{
public:
    ExtForceField() : fpField_(nullptr) {}

    // From IMDModule
    IMdpOptionProvider* mdpOptionProvider() override { return this; }
    IMDOutputProvider*  outputProvider() override { return this; }
    void                initForceProviders(ForceProviders* forceProviders) override
    {
        if (isActive())
        {
            forceProviders->addForceProvider(this);
        }
    }

    // From IMdpOptionProvider
    void initMdpTransform(IKeyValueTreeTransformRules* transform) override;
    void initMdpOptions(IOptionsContainerWithSections* options) override;
    void buildMdpOutput(KeyValueTreeObjectBuilder* builder) const override;

    // From IMDOutputProvider
    void initOutput(FILE* fplog, int nfile, const t_filenm fnm[], bool bAppendFiles, const gmx_output_env_t* oenv) override;
    void finishOutput() override;

    // From IForceProvider
    //! \copydoc IForceProvider::calculateForces()
    void calculateForces(const ForceProviderInput& forceProviderInput,
                         ForceProviderOutput*      forceProviderOutput) override;

    void subscribeToSimulationSetupNotifications(MdModulesNotifier* /* notifier */) override {}
    void subscribeToPreProcessingNotifications(MdModulesNotifier* /* notifier */) override {}

private:
    //! Return whether or not to apply a field
    bool isActive() const;

    /*! \brief Return the field strength
     *
     * \param[in] dim The spatial direction
     * \param[in] t   The time (ps)
     * \return The field strength in kJ/mol/nm units
     */
    real field(int dim, real t) const;

    /*! \brief Print the field components to a file
     *
     * \param[in] t   The time
     * Will throw and exit with fatal error if file is not open.
     */
    void printComponents(double t) const;

    //! The components of the applied extforce field in each coordinate dimension
    ExtForceFieldDimension efield_[DIM];
    //! File pointer for extforce field
    FILE* fpField_;
};

//! Converts dynamic parameters from new mdp format to (E0, omega, t0, sigma, zmin, zmax).
void convertParameters(gmx::KeyValueTreeObjectBuilder* builder, const std::string& value)
{
    const std::vector<std::string> sxt = splitString(value);
    if (sxt.empty())
    {
        return;
    }
    if (sxt.size() != 6)
    {
        GMX_THROW(InvalidInputError("Please specify E0 omega t0 sigma zmin zmax for extforce fields"));
    }
    builder->addValue<real>("E0", fromString<real>(sxt[0]));
    builder->addValue<real>("omega", fromString<real>(sxt[1]));
    builder->addValue<real>("t0", fromString<real>(sxt[2]));
    builder->addValue<real>("sigma", fromString<real>(sxt[3]));
    builder->addValue<real>("zmin", fromString<real>(sxt[4]));
    builder->addValue<real>("zmax", fromString<real>(sxt[5]));
}

void ExtForceField::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    rules->addRule().from<std::string>("/extforce-field-x").toObject("/extforce-field/x").transformWith(&convertParameters);
    rules->addRule().from<std::string>("/extforce-field-y").toObject("/extforce-field/y").transformWith(&convertParameters);
    rules->addRule().from<std::string>("/extforce-field-z").toObject("/extforce-field/z").transformWith(&convertParameters);
}

void ExtForceField::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection("extforce-field"));
    efield_[XX].initMdpOptions(&section, "x");
    efield_[YY].initMdpOptions(&section, "y");
    efield_[ZZ].initMdpOptions(&section, "z");
}

void ExtForceField::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    std::string comment =
            R"(; ExtForce fields
; Format for extforce-field-x, etc. is: six real variables:
; amplitude (V/nm), frequency omega (1/ps), time for the pulse peak (ps),
; width of the pulse sigma (ps), the lower z (nm) and highest z (nm) to consider.
; Omega = 0 means static field,
; sigma = 0 means no pulse, leaving the field to be a cosine function,
; zmin=zmax means that the field is applied in whole system.)";
    builder->addValue<std::string>("comment-extforce-field", comment);
    efield_[XX].buildMdpOutput(builder, "x");
    efield_[YY].buildMdpOutput(builder, "y");
    efield_[ZZ].buildMdpOutput(builder, "z");
}

void ExtForceField::initOutput(FILE* fplog, int nfile, const t_filenm fnm[], bool bAppendFiles, const gmx_output_env_t* oenv)
{
    if (isActive())
    {
        please_cite(fplog, "Caleman2008a");

        // Optional output file showing the field, see manual.
        if (opt2bSet("-field", nfile, fnm))
        {
            if (bAppendFiles)
            {
                fpField_ = gmx_fio_fopen(opt2fn("-field", nfile, fnm), "a+");
            }
            else
            {
                fpField_ = xvgropen(opt2fn("-field", nfile, fnm), "Applied extforce field",
                                    "Time (ps)", "E (kJ/mol/nm)", oenv);
            }
        }
    }
}

void ExtForceField::finishOutput()
{
    if (fpField_ != nullptr)
    {
        // This is opened sometimes with xvgropen, sometimes with
        // gmx_fio_fopen, so we use the least common denominator for closing.
        gmx_fio_fclose(fpField_);
        fpField_ = nullptr;
    }
}

real ExtForceField::field(int dim, real t) const
{
    return efield_[dim].evaluate(t);
}


bool ExtForceField::isActive() const
{
    return (efield_[XX].a() != 0 || efield_[YY].a() != 0 || efield_[ZZ].a() != 0);
}

void ExtForceField::printComponents(double t) const
{
    fprintf(fpField_, "%10g  %10g  %10g  %10g\n", t, field(XX, t), field(YY, t), field(ZZ, t));
}

void ExtForceField::calculateForces(const ForceProviderInput& forceProviderInput,
                                    ForceProviderOutput*      forceProviderOutput)
{
    // GROUP works like this :
    // for "user1-grps              = NA CL OW"
    // NA -> 0
    // CL -> 1
    // OW -> 2
    // All others -> 3
    if (isActive())
    {
        const t_mdatoms& mdatoms = forceProviderInput.mdatoms_;
        const ArrayRef<const RVec> x = forceProviderInput.x_;
        const double     t       = forceProviderInput.t_;
        const t_commrec& cr      = forceProviderInput.cr_;

        // NOTE: The non-conservative extforce field does not have a virial
        ArrayRef<RVec> f = forceProviderOutput->forceWithVirial_.force_;

        for (int m = 0; (m < DIM); m++)
        {   
            const real zmin = efield_[m].zmin();
            const real zmax = efield_[m].zmax();
            const real fieldStrength = field(m, t);

            if (fieldStrength == 0)
            {
                continue;
            }
            
            // Sum up total mass for atoms such as zmin<=z<zmax
            real mtot1 = 0; 
            for (int i = 0; i < mdatoms.homenr; ++i)
            {
                // Skip atoms that do not belong to USER1 group (if one exists)
                if (mdatoms.cU1!=NULL){
                    if (mdatoms.cU1[i]!=0){
                        continue;
                    }
                }
                // If zmin<zmax, include atoms that belong to this defined area
                if (x[i][2] >= zmin && x[i][2] < zmax)
                {
                    mtot1 += mdatoms.massT[i];
                }
                // If zmin==zmax, the force is applied on atoms without condition on their location 
                else if (zmin == zmax)
                {
                    mtot1 += mdatoms.massT[i];
                }
            }

            // Assign forces to atoms
            if (mtot1 != 0){
                for (int i = 0; i < mdatoms.homenr; ++i)
                {
                    // Skip atoms that do not belong to USER1 group (if one exists)
                    if (mdatoms.cU1!=NULL){
                        if (mdatoms.cU1[i]!=0){
                            continue;
                        }
                    }
                    // If zmin<zmax, include atoms that belong to this defined area
                    if (x[i][2] >= zmin && x[i][2] < zmax)
                    {
                        f[i][m] += mdatoms.massT[i]*fieldStrength/mtot1;
                    }
                    // If zmin==zmax, the force is applied on atoms without condition on their location 
                    else if (zmin == zmax)
                    {
                        f[i][m] += mdatoms.massT[i]*fieldStrength/mtot1;
                    }
                }
            }

            if (mdatoms.cU2!=NULL){
                // Sum up total mass for atoms such as zmin<=z<zmax
                real mtot2 = 0; 
                for (int i = 0; i < mdatoms.homenr; ++i)
                {
                    // Skip atoms that do not belong to USER2 group
                    if (mdatoms.cU2[i]!=0){
                        continue;
                    }
                    // If zmin<zmax, include atoms that belong to this defined area
                    if (x[i][2] >= zmin && x[i][2] < zmax)
                    {
                        mtot2 += mdatoms.massT[i];
                    }
                    // If zmin==zmax, the force is applied on atoms without condition on their location 
                    else if (zmin == zmax)
                    {
                        mtot2 += mdatoms.massT[i];
                    }
                }

                // Assign forces to atoms
                if (mtot2 != 0){
                    for (int i = 0; i < mdatoms.homenr; ++i)
                    {
                        // Skip atoms that do not belong to USER2 group
                        if (mdatoms.cU2[i]!=0){
                            continue;
                        }
                        // If zmin<zmax, include atoms that belong to this defined area
                        if (x[i][2] >= zmin && x[i][2] < zmax)
                        {
                            f[i][m] += -mdatoms.massT[i]*fieldStrength/mtot2;
                        }
                        // If zmin==zmax, the force is applied on atoms without condition on their location 
                        else if (zmin == zmax)
                        {
                            f[i][m] += -mdatoms.massT[i]*fieldStrength/mtot2;
                        }
                    }
                }
            }

        }
        if (MASTER(&cr) && fpField_ != nullptr)
        {
        }
    }
}

} // namespace

std::unique_ptr<IMDModule> createExtForceFieldModule()
{
    return std::make_unique<ExtForceField>();
}

} // namespace gmx
