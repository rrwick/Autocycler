#!/bin/bash

install_tool() {
    set -e # Exit immediately if a command exits with a non-zero status.

    if [ "$#" -lt 1 ]; then
        echo "Error: Usage: $0 <env_name> [package1 package2 ...]" >&2
        echo "Error: At least the environment name is required." >&2
        return 1
    fi

    local env_name="$1"
    shift

    local primary_package_to_install
    local additional_packages_to_install=()
    local packages_for_echo_list=()

    if [ "$#" -eq 0 ]; then
        primary_package_to_install="$env_name"
        packages_for_echo_list=("$env_name")
    else
        primary_package_to_install="$1"
        packages_for_echo_list+=("$primary_package_to_install")
        shift
        additional_packages_to_install=("$@")

        if [ "${#additional_packages_to_install[@]}" -gt 0 ]; then
            packages_for_echo_list+=("${additional_packages_to_install[@]}")
        fi
    fi

    echo ">>> Creating environment '$env_name' and installing: ${packages_for_echo_list[*]}"
    micromamba create -y -n "$env_name" -c conda-forge -c bioconda "$primary_package_to_install" "${additional_packages_to_install[@]}"

    echo ">>> Successfully created environment '$env_name' and installed ${packages_for_echo_list[*]}"
    micromamba clean --all --yes
}

# Call the function with all script arguments
install_tool "$@"
