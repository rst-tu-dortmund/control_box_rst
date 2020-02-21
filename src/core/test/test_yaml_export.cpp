/*********************************************************************
 *
 *  Software License Agreement
 *
 *  Copyright (c) 2020,
 *  TU Dortmund - Institute of Control Theory and Systems Engineering.
 *  All rights reserved.
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 *  Authors: Christoph Rösmann
 *********************************************************************/

#ifdef YAML_SUPPORT

#include <corbo-core/macros.h>

#include <corbo-core/signals.h>
#include <corbo-core/yaml_export.h>

#include "gtest/gtest.h"



class TestYAMLExport : public testing::Test
{
 protected:
    // You can do set-up work for each test here.
    TestYAMLExport() {}
    // You can do clean-up work that doesn't throw exceptions here.
    virtual ~TestYAMLExport() {}
    // If the constructor and destructor are not enough for setting up
    // and cleaning up each test, you can define the following methods:

    // Code here will be called immediately after the constructor (right
    // before each test).
    // virtual void SetUp() {}
    // Code here will be called immediately after each test (right
    // before the destructor).
    // virtual void TearDown();
};

TEST_F(TestYAMLExport, test1)
{
    corbo::YamlExporter yaml_exporter;
    corbo::TimeSeriesSignal ts_signal;

    std::string filename = "testfile.yaml";

    EXPECT_TRUE(yaml_exporter.exportTimeSeriesSignal(filename, ts_signal));
}

TEST_F(TestYAMLExport, export_matrix)
{
    corbo::YamlExporter yaml_exporter;
    std::string filename = "testfile.yaml";

    Eigen::Matrix3d mat = Eigen::Matrix3d::Identity();
    mat(0, 1) = 2;  // for testing column/row major
    corbo::MatrixSignal mat_sig(mat, "test_matrix");
    mat_sig.header.name = "test_header_name";
    mat_sig.header.time = corbo::Time(1);

    EXPECT_TRUE(yaml_exporter.exportMatrixSignal(filename, mat_sig));
}

TEST_F(TestYAMLExport, export_matrix_set)
{
    corbo::YamlExporter yaml_exporter;
    std::string filename = "testfile.yaml";

    Eigen::Matrix3d mat = Eigen::Matrix3d::Identity();
    mat(0, 1) = 2;  // for testing column/row major
    corbo::MatrixSignal::Ptr mat_sig = std::make_shared<corbo::MatrixSignal>(mat, "test_matrix");
    mat_sig->header.name             = "mat_sig_name";
    mat_sig->header.time             = corbo::Time(1);

    corbo::MatrixSetSignal mat_set_sig;
    mat_set_sig.header.name = "mat_set_sig_name";
    mat_sig->header.time    = corbo::Time(2);

    mat_set_sig.add(mat_sig);
    mat_set_sig.add(mat_sig);

    Eigen::Matrix2d mat2d;
    mat2d << 1, 2, 3, 4;
    mat_set_sig.add(mat2d, "2d_mat");

    EXPECT_TRUE(yaml_exporter.exportMatrixSetSignal(filename, mat_set_sig));
}

#endif // YAML_SUPPORT
