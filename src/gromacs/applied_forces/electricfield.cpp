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
 * Defines data structure and utilities for electric fields
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include "electricfield.h"

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
 * \brief Describes an applied electric field in a coordinate
 * dimension.
 *
 * Can compute the applied field strength at a time, and supports
 * operations to read and form corresponding .mdp contents.
 */
class ElectricFieldDimension
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
        builder->addUniformArray<real>("electric-field-" + name, { a_, omega_, t0_, sigma_, zmin_, zmax_ });
    }

    /*! \brief Evaluates this field component at given time.
     *
     * \param[in] t The time to evualate at
     * \return The electric field
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

    /*! \brief Evaluates this field component at given time and z coordinate.
     *
     * \param[in] t The time to evualate at
     * \param[in] z z coordinate
     * \return The electric field
     */
    real evaluate_at_z(real t, real z) const
    {
        if (zmin_ == zmax_)
        {
            return 0;
        }
        else
        {
            return evaluate(t);
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
    //! Higher value of z coordinate (nm) for which the field is applied.
    real zmax_ = 0;
};

/*! \internal
 * \brief Describe time dependent electric field
 *
 * Class that implements a force to be evaluated in mdrun.
 * The electric field can be pulsed and oscillating, simply
 * oscillating, or static, in each of X,Y,Z directions.
 */
class ElectricField final : public IMDModule, public IMdpOptionProvider, public IMDOutputProvider, public IForceProvider
{
public:
    ElectricField() : fpField_(nullptr) {}

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
     * \return The field strength in V/nm units
     */
    real field(int dim, real t) const;

    /*! \brief Return the field strength
     *
     * \param[in] dim The spatial direction
     * \param[in] t   The time (ps)
     * \param[in] z   The position along z
     * \return The field strength in V/nm units
     */
    real field_at_z(int dim, real t, real z) const;

    /*! \brief Print the field components to a file
     *
     * \param[in] t   The time
     * Will throw and exit with fatal error if file is not open.
     */
    void printComponents(double t) const;

    //! The components of the applied electric field in each coordinate dimension
    ElectricFieldDimension efield_[DIM];
    //! File pointer for electric field
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
        GMX_THROW(InvalidInputError("Please specify E0 omega t0 sigma zmin zmax for electric fields"));
    }
    builder->addValue<real>("E0", fromString<real>(sxt[0]));
    builder->addValue<real>("omega", fromString<real>(sxt[1]));
    builder->addValue<real>("t0", fromString<real>(sxt[2]));
    builder->addValue<real>("sigma", fromString<real>(sxt[3]));
    builder->addValue<real>("zmin", fromString<real>(sxt[4]));
    builder->addValue<real>("zmax", fromString<real>(sxt[5]));
}

void ElectricField::initMdpTransform(IKeyValueTreeTransformRules* rules)
{
    rules->addRule().from<std::string>("/electric-field-x").toObject("/electric-field/x").transformWith(&convertParameters);
    rules->addRule().from<std::string>("/electric-field-y").toObject("/electric-field/y").transformWith(&convertParameters);
    rules->addRule().from<std::string>("/electric-field-z").toObject("/electric-field/z").transformWith(&convertParameters);
}

void ElectricField::initMdpOptions(IOptionsContainerWithSections* options)
{
    auto section = options->addSection(OptionSection("electric-field"));
    efield_[XX].initMdpOptions(&section, "x");
    efield_[YY].initMdpOptions(&section, "y");
    efield_[ZZ].initMdpOptions(&section, "z");
}

void ElectricField::buildMdpOutput(KeyValueTreeObjectBuilder* builder) const
{
    std::string comment =
            R"(; Electric fields
; Format for electric-field-x, etc. is: four real variables:
; amplitude (V/nm), frequency omega (1/ps), time for the pulse peak (ps),
; width of the pulse sigma (ps), the lower z (nm) and highest z (nm) to consider.
; Omega = 0 means static field,
; sigma = 0 means no pulse, leaving the field to be a cosine function,
; zmin=zmax means that the field is applied in whole system.)";
    efield_[XX].buildMdpOutput(builder, "x");
    efield_[YY].buildMdpOutput(builder, "y");
    efield_[ZZ].buildMdpOutput(builder, "z");
}

void ElectricField::initOutput(FILE* fplog, int nfile, const t_filenm fnm[], bool bAppendFiles, const gmx_output_env_t* oenv)
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
                fpField_ = xvgropen(opt2fn("-field", nfile, fnm), "Applied electric field",
                                    "Time (ps)", "E (V/nm)", oenv);
            }
        }
    }
}

void ElectricField::finishOutput()
{
    if (fpField_ != nullptr)
    {
        // This is opened sometimes with xvgropen, sometimes with
        // gmx_fio_fopen, so we use the least common denominator for closing.
        gmx_fio_fclose(fpField_);
        fpField_ = nullptr;
    }
}

real ElectricField::field(int dim, real t) const
{
    return efield_[dim].evaluate(t);
}

real ElectricField::field_at_z(int dim, real t, real z) const
{
    return efield_[dim].evaluate_at_z(t, z);
}

bool ElectricField::isActive() const
{
    return (efield_[XX].a() != 0 || efield_[YY].a() != 0 || efield_[ZZ].a() != 0);
}

void ElectricField::printComponents(double t) const
{
    fprintf(fpField_, "%10g  %10g  %10g  %10g\n", t, field(XX, t), field(YY, t), field(ZZ, t));
}

void ElectricField::calculateForces(const ForceProviderInput& forceProviderInput,
                                    ForceProviderOutput*      forceProviderOutput)
{
    if (isActive())
    {
        const t_mdatoms& mdatoms = forceProviderInput.mdatoms_;
        const ArrayRef<const RVec> x = forceProviderInput.x_;
        const double     t       = forceProviderInput.t_;
        const t_commrec& cr      = forceProviderInput.cr_;

        // NOTE: The non-conservative electric field does not have a virial
        ArrayRef<RVec> f = forceProviderOutput->forceWithVirial_.force_;

        for (int m = 0; (m < DIM); m++)
        {
            const real zmin = efield_[m].zmin();
            const real zmax = efield_[m].zmax();
            const real fieldStrength = FIELDFAC * field(m, t);
            
            if (fieldStrength == 0)
            {
                continue;
            }

            for (int i = 0; i < mdatoms.homenr; ++i)
            {
                if (x[i][2] >= zmin && x[i][2] < zmax)
                {
                    f[i][m] += mdatoms.chargeA[i] * fieldStrength;
                }
                else if (zmin == zmax)
                {
                    f[i][m] += mdatoms.chargeA[i] * fieldStrength;
                }
            }
        }
        if (MASTER(&cr) && fpField_ != nullptr)
        {
            printComponents(t);
        }
    }
}

} // namespace

std::unique_ptr<IMDModule> createElectricFieldModule()
{
    return std::make_unique<ElectricField>();
}

} // namespace gmx
