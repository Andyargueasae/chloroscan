# !/bin/bash

# Check if FragGeneScanRs is installed and added to the PATH.

check_fraggenescans() {
    if command -v FragGeneScanRs &> /dev/null
    then
        echo "FragGeneScanRs is already installed and added to the PATH."
        return 0
    else
        echo "FragGeneScanRs is not installed or not in the PATH."
        return 1 # Forward to installation of FragGeneScanRs.
    fi
}

install_fraggenescans() {
    # Firstly check if the cargo is in command.
    if command -v cargo &> /dev/null
    then
        echo "Cargo is already installed."
    else
        echo "Cargo is not installed."
        # The program assumes that the user has the curl installed.
        curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
        source $HOME/.cargo/env
    fi

    cargo install frag_gene_scan_rs

    CARGO_BIN_DIR="$HOME/.cargo/bin"

    if [[ ":$PATH:" != *":$CARGO_BIN_DIR:"* ]]; then
        echo "Adding Cargo bin directory ($CARGO_BIN_DIR) to PATH."
        export PATH=$CARGO_BIN_DIR:$PATH

        # Persist this PATH update in .bashrc
        if ! grep -q "$CARGO_BIN_DIR" ~/.bashrc; then
            echo "export PATH=$CARGO_BIN_DIR:\$PATH" >> ~/.bashrc
            echo "Added Cargo bin directory to PATH in ~/.bashrc."
        fi
    fi

    echo "FragGeneScanRs has been installed and added to the PATH."
    # Add fraggenescans to the PATH.

}

if ! check_fraggenescans; then
    install_fraggenescans
else
    echo "No installation of FragGeneScanRs required."
fi





