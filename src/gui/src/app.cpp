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

#include <corbo-gui/main_window.h>

#include <QApplication>

#if defined(_WIN32) && defined(_MINGW)
#include <QtPlugin>
Q_IMPORT_PLUGIN(QWindowsIntegrationPlugin)
#endif  // defined

void parseArguments(int argc, char** argv, bool* show_console);

int main(int argc, char** argv)
{
    bool show_console = false;
    parseArguments(argc, argv, &show_console);

    QApplication app(argc, argv);
    // set organization and application name which is required for storing app data
    app.setOrganizationName(corbo::gui::OrganizationName);
    app.setApplicationName(corbo::gui::ApplicationName);

    corbo::gui::corboMainWindow main_window;
    main_window.show();

#ifdef _WIN32
    if (!show_console) FreeConsole();
#endif

    // enter main application event-loop
    return app.exec();
}

void parseArguments(int argc, char** argv, bool* show_console)
{
    if (show_console) *show_console = false;

    for (int i = 1; i < argc; ++i)
    {
        if (strcmp(argv[i], "-c") == 0 || strcmp(argv[i], "--c") == 0 || strcmp(argv[i], "-console") == 0)  // show console
        {
            if (show_console) *show_console = true;
        }
    }
}
