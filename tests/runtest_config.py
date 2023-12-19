def configure(options, input_files,extra_arg):
    """
    This function is used by runtest to configure runtest
    at runtime for code specific launch command and file naming.
    """

    from os import path
    from sys import platform

    launcher = 'FRET_Embedlab'
    launcher_full_path = path.normpath(path.join(options.binary_dir, launcher))

    if len(input_files) == 1:
        (inp) = input_files
        extra = None

    inp_no_prefix = inp[0]
    command = []
    command.append(launcher_full_path)
    command.append(inp_no_prefix)

    full_command = ' '.join(command)
    print(full_command)

    output_prefix = inp_no_prefix[:-4]
    print(output_prefix)

    relative_reference_path = 'reference'

    return launcher, full_command, output_prefix, relative_reference_path
