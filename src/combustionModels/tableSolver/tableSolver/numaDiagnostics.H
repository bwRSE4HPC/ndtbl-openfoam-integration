// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef NDTBL_NUMA_DIAGNOSTICS_H_
#define NDTBL_NUMA_DIAGNOSTICS_H_

#include "ndtbl/diagnostics.hpp"

#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Foam
{
namespace combustionModels
{
namespace detail
{

inline const char* numaTsvHeader()
{
    return "# ndtbl_numa_tsv_version=1\n"
        "sample\ttime\tlabel\tgroup\trank\thost\tnuma_available\tpolicy"
        "\tkernel_page_size_bytes\tnuma_node\tpages\n";
}

/// One row per reported NUMA node, or a placeholder for missing node counters.
inline std::string numaTsvRows
(
    const std::size_t sample,
    const double time,
    const std::string& label,
    const std::size_t group,
    const int rank,
    const std::string& host,
    const ndtbl::residency_info& info
)
{
    for (const auto* value : {&label, &host, &info.numa_maps_policy})
    {
        if (value->find_first_of("\t\r\n") != std::string::npos)
        {
            throw std::invalid_argument("NUMA TSV field contains a tab or newline");
        }
    }

    std::ostringstream prefix;
    prefix << std::setprecision(std::numeric_limits<double>::max_digits10)
        << sample << '\t' << time << '\t' << label << '\t' << group
        << '\t' << rank << '\t' << host
        << '\t' << (info.numa_maps_available ? 1 : 0)
        << '\t' << (info.numa_maps_available ? info.numa_maps_policy : "nan")
        << '\t';
    if (info.numa_maps_available && info.numa_maps_kernel_page_size_available)
    {
        prefix << info.numa_maps_kernel_page_size_bytes;
    }
    else
    {
        prefix << "nan";
    }

    if (!info.numa_maps_available || info.numa_maps_node_pages.empty())
    {
        return prefix.str() + "\tnan\tnan\n";
    }

    std::ostringstream rows;
    for (const auto& node : info.numa_maps_node_pages)
    {
        rows << prefix.str() << '\t' << node.first << '\t' << node.second << '\n';
    }
    return rows.str();
}

} // namespace detail
} // namespace combustionModels
} // namespace Foam

#endif
