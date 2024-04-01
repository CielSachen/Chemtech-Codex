from functools import reduce
from math import fsum, log10
from os import path, remove
from typing import Self

from prettytable import PrettyTable

# Defines the path to where the file containing the results would be created.
file_path = "results.txt"


def subscript(formula: str) -> str:
    """Converts all the integers in a string of a chemical formula into its subscript form.

    Args:
        formula (str): The string of a formula.

    Returns:
        str: The string of a chemical formula with subscripts.
    """
    # Create a translation dictionary for the integers and their subscript form.
    subscript_translation = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

    # Separates the parts of the given chemical formula into its coefficient and main.
    coefficient = formula[0]
    main_formula = formula[1:]

    # Translates the integers in the main chemical formula into their subscript forms.
    main_formula = main_formula.translate(subscript_translation)

    # Returns the string of a chemical formula with subscripts.
    return coefficient + main_formula


def superscript(integers: str) -> str:
    """Converts all the integers in a string into its subscript form.

    Args:
        formula (str): The string.

    Returns:
        str: The string with subscripts.
    """
    # Create a translation dictionary for the integers and their superscript form.
    superscript_translation = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")

    # Translates the integers into their superscript forms and returns it.
    return integers.translate(superscript_translation)


class Thermochemistry:
    # Creates a new ASCII table to store the givens.
    givens_table = PrettyTable(
        ["Formula", "Type", "Number", "Standard Heat (ΔH°f)", "Standard Entropy (S°)"]
    )

    def __init__(self: Self, raw_chemical_equation: str) -> None:
        """Constructs a new instance of the Thermochemistry class,

        Args:
            self (Self): The class of the calculator.
            raw_chemical_equation (str): The chemical equation.
        """
        # Separates the chemical equation into its reactants and products.
        chemical_equation = raw_chemical_equation.replace(" ", "").split("->")
        # Formats and splits the reactants and products into individual formulas.
        self.reactants = list(map(subscript, chemical_equation[0].split("+")))
        self.products = list(map(subscript, chemical_equation[1].split("+")))
        # Rejoins the reactants and products into a stylized chemical equation.
        self.chemical_equation = " → ".join(
            [" + ".join(self.reactants), " + ".join(self.products)]
        )

    def get_thermodynamic_values(
        self: Self, formulas: list[str], formula_type: str
    ) -> tuple[list[float], list[float]]:
        """Gets the thermodynamic values of chemical formulas through inputs from the user.

        Args:
            self (Self): The class of the calculator.
            formulas (list[str]): The chemical formulas.
            formula_type (str): The type of the thermodynamic values.

        Returns:
            tuple[list[float], list[float]]: The user given standard heats and entropies.
        """
        # Create lists to store the standard heats and entropies of the given formulas.
        standard_heats: list[float] = []
        standard_entropies: list[float] = []

        # Loop through each of the given chemical formulas.
        for formula in formulas:
            # Resets the variable to allow for looping if inputted value is invalid
            successful = False

            # Loops the input prompt until the user gives a valid thermodynamic value.
            while not successful:
                # Tries to prompt the user to input the thermodynamic values.
                try:
                    # Checks if a chemical formula has a coefficient.
                    if formula[0].isdigit():
                        # Separates the parts of the given chemical formula into its coefficient and main.
                        coefficient = int(formula[0])
                        main_formula = formula[1:]

                        # Asks the user for the standard heat value and multiplies it to the coefficient.
                        given_standard_heat = float(
                            input(
                                f'\nPlease input the STANDARD HEAT (ΔH°f) of "{main_formula}" (float or int ONLY): '
                            )
                        )
                        standard_heat = given_standard_heat * coefficient
                        # Asks the user for the standard entropy value and multiplies it to the coefficient.
                        given_standard_entropy = float(
                            input(
                                f'Please input the STANDARD ENTROPY (S°) of "{main_formula}" (float or int ONLY): '
                            )
                        )
                        standard_entropy = given_standard_entropy * coefficient

                        # Adds the given thermodynamic values to the table of givens
                        self.givens_table.add_row(  # type: ignore
                            [
                                main_formula,
                                formula_type,
                                coefficient,
                                f"{given_standard_heat} kJ/mol",
                                f"{given_standard_entropy} J/K × mol",
                            ]
                        )
                    else:
                        # Asks the user for the standard heat value.
                        standard_heat = given_standard_heat = float(
                            input(
                                f'\nPlease input the STANDARD HEAT (ΔH°f) of "{formula}" (float or int ONLY): '
                            )
                        )
                        # Asks the user for the standard entropy value.
                        standard_entropy = given_standard_entropy = float(
                            input(
                                f'Please input the STANDARD ENTROPY (S°) of "{formula}" (float or int ONLY): '
                            )
                        )

                        # Adds the given thermodynamic values to the table of givens
                        self.givens_table.add_row(  # type: ignore
                            [
                                formula,
                                formula_type,
                                1,
                                f"{given_standard_heat} kJ/mol",
                                f"{given_standard_entropy} J/K × mol",
                            ]
                        )

                    # Adds the given standard heat and entropy to the lists that will be returned by the function
                    standard_heats.append(standard_heat)
                    standard_entropies.append(standard_entropy)

                    # Stops the input prompt from repeating.
                    successful = True
                # Catches error exception for when the given value is invalid.
                except ValueError:
                    # Informs the user what values are valid.
                    print("\nThe given value must be a valid float or int ONLY.")

                    # Repeats the input prompt.
                    successful = False

        # Returns the lists of given standard heat and entropy values.
        return standard_heats, standard_entropies

    def calculate_rxn_enthalpy(
        self: Self,
        product_standard_heats: list[float],
        reactant_standard_heats: list[float],
    ) -> float:
        """Calculates the standard enthalpy of the standard heats of the products and reactants.

        Args:
            self (Self): The class of the calculator.
            product_standard_heats (list[float]): The standard heats of the products.
            reactant_standard_heats (list[float]): The standard heats of the reactants.

        Returns:
            float: The calculated standard enthalpy.
        """
        # Opens the file where the results are contained.
        file = open(file_path, "a", encoding="utf-8")

        # Writes to the file the equation that will be used to calculate the standard enthalpy.
        file.write(
            "\n\nEquation for standard enthalpy (ΔH°rxn):\nΔH°rxn = Σ(nΔS° products) - Σ(nΔS° reactants)"
        )

        # Writes to the file the solution for calculating the standard enthalpy.
        file.write(
            f'\n\nSolution for the standard enthalpy (ΔH°rxn) of "{self.chemical_equation}":'
        )
        file.write(
            f"\nΔH°rxn = ({' + '.join(map(str, product_standard_heats))}) - ({' + '.join(map(str, reactant_standard_heats))}) kJ/mol"
        )

        # Calculates the standard enthalpy.
        standard_enthalpy = fsum(product_standard_heats) - fsum(reactant_standard_heats)

        # Writes to the file the calculated standard enthalpy.
        file.write(f"\nΔH°rxn = {standard_enthalpy} kJ/mol")

        # Checks if the enthalpy is greater than zero then writes to the file that the reaction is endothermic.
        if standard_enthalpy > 0:
            file.write(
                "\n\nBecause the standard enthalpy (ΔH°rxn) is greater than 0, the reaction must be ENDOTHERMIC."
            )
        # Checks if the enthalpy is less than zero then writes to the file that the reaction is exothermic.
        elif standard_enthalpy < 0:
            file.write(
                "\n\nBecause the standard enthalpy (ΔH°rxn) is less than 0, the reaction must be EXOTHERMIC."
            )
        # Writes to the file that the reaction is neither endothermic or exothermic.
        else:
            file.write(
                "\n\nBecause the standard enthalpy (ΔH°rxn) is equal to 0, the reaction is neither endothermic or exothermic."
            )

        # Closes the opened file.
        file.close()

        # Returns the calculated standard enthalpy.
        return standard_enthalpy

    def calculate_rxn_entropy(
        self: Self,
        product_standard_entropies: list[float],
        reactant_standard_entropies: list[float],
    ) -> float:
        """Calculates the standard entropy change of the standard entropies of the products and reactants.

        Args:
            self (Self): The class of the calculator.
            product_standard_entropies (list[float]): The standard entropies of the products.
            reactant_standard_entropies (list[float]): The standard entropies of the reactants.

        Returns:
            float: The calculated standard entropy change.
        """
        # Opens the file where the results are contained.
        file = open(file_path, "a", encoding="utf-8")

        # Writes to the file the equation that will be used to calculate the standard entropy change.
        file.write(
            "\n\nEquation for standard entropy change / entropy of the system (ΔS°rxn / ΔS°sys):\nΔS°rxn = Σ(nΔH°f products) - Σ(nΔH°f reactants)"
        )

        # Writes to the file the solution for calculating the standard entropy change.
        file.write(
            f'\n\nSolution for the standard entropy change / entropy of the system (ΔS°rxn / ΔS°sys) of "{self.chemical_equation}":'
        )
        file.write(
            f"\nΔS°rxn = {' + '.join(map(str, product_standard_entropies))} - {' + '.join(map(str, reactant_standard_entropies))} J/K × mol"
        )

        # Calculates the standard entropy change.
        standard_entropy_change = fsum(product_standard_entropies) - fsum(
            reactant_standard_entropies
        )

        # Writes to the file the calculated standard entropy change.
        file.write(f"\nΔS°rxn = {standard_entropy_change} J/K × mol")

        # Closes the opened file.
        file.close()

        # Returns the calculated standard entropy change.
        return standard_entropy_change

    def calculate_surr_entropy(
        self: Self, system_entropy: float, temperature: float
    ) -> float:
        """Calculates the enthalpy of the surroundings of the system entropy and temperature.

        Args:
            self (Self): The class of the calculator.
            system_entropy (list[float]): The system entropy or standard entropy change.
            temperature (list[float]): The temperature.

        Returns:
            float: The calculated entropy of the surroundings.
        """
        # Opens the file where the results are contained.
        file = open(file_path, "a", encoding="utf-8")

        # Writes to the file the equation that will be used to calculate the entropy of the surroundings.
        file.write(
            "\n\nEquation for entropy of the surroundings (ΔS°surr):\nΔS°surr = -(ΔS°sys) / T"
        )

        # Writes to the file the solution for calculating the entropy of the surroundings.
        file.write(
            f'\n\nSolution for the entropy of the surroundings (ΔS°surr) of "{self.chemical_equation} at {temperature} K":'
        )
        file.write(f"\nΔS°surr = -({system_entropy}) J/mol / {temperature} K")

        # Calculates the entropy of the surroundings.
        surrounding_entropy = (system_entropy * -1) / temperature

        # Writes to the file the calculated entropy of the surroundings.
        file.write(f"\nΔS°surr = {surrounding_entropy} J/K × mol")

        # Closes the opened file.
        file.close()

        # Returns the calculated entropy of the surroundings.
        return surrounding_entropy

    def calculate_univ_entropy(
        self: Self,
        system_entropy: float,
        surrounding_entropy: float,
    ) -> float:
        """Calculates the entropy of the universe of the the system entropy and entropy of the surroundings.

        Args:
            self (Self): The class of the calculator.
            system_entropy (list[float]): The system entropy or standard entropy change.
            surrounding_entropy (list[float]): The entropy of the surroundings.

        Returns:
            float: The calculated entropy of the universe.
        """
        # Opens the file where the results are contained.
        file = open(file_path, "a", encoding="utf-8")

        # Writes to the file the equation that will be used to calculate the entropy of the universe.
        file.write(
            "\n\nEquation for entropy of the universe (ΔS°univ):\nΔS°univ = ΔS°sys ΔS°surr"
        )

        # Writes to the file the solution for calculating the entropy of the universe.
        file.write(
            f'\n\nSolution for the entropy of the universe (ΔS°univ) of "{self.chemical_equation}":'
        )
        file.write(
            f"\nΔS°univ = {system_entropy} J/K × mol + {surrounding_entropy} J/K × mol"
        )

        # Calculates the entropy of the universe.
        universal_entropy = system_entropy + surrounding_entropy

        # Writes to the file the calculated entropy of the universe.
        file.write(f"\nΔS°univ = {universal_entropy} J/K × mol")

        # Checks if the enthalpy is greater than zero then writes to the file that the reaction is spontaneous.
        if universal_entropy > 0:
            file.write(
                "\n\nBecause the entropy of the universe (ΔS°univ) is greater than 0, the reaction must be SPONTANEOUS."
            )
        # Checks if the enthalpy is less than zero then writes to the file that the reaction is non-spontaneous.
        elif universal_entropy < 0:
            file.write(
                "\n\nBecause the entropy of the universe (ΔS°univ) is less than 0, the reaction must be NON-SPONTANEOUS."
            )
        # Writes to the file that the reaction is at equilibrium.
        else:
            file.write(
                "\n\nBecause the entropy of the universe (ΔS°univ) is equal to 0, the reaction is at equilibrium."
            )

        # Closes the opened file.
        file.close()

        # Returns the calculated entropy of the universe.
        return universal_entropy

    def calculate(self: Self) -> None:
        """Runs the thermodynamics calculator.

        Args:
            self (Self): The class of the calculator.
        """
        # Declares the default temperature value.
        temperature = 298.0

        # Resets the variable to allow for looping if inputted value is invalid
        successful = False

        # Loops the input prompt until the user gives a valid temperature value.
        while not successful:
            # Tries to prompt the user to input the temperature value.
            try:
                # Asks the user for the temperature that the chemical reaction occurs in.
                given_temperature = input(
                    "Please input the temperature in Kelvins (float, int, or [press ENTER for 298] ONLY): "
                )

                # Checks if the user gave a temperature
                if not given_temperature == "":
                    # Sets the temperature to the given temperature.
                    temperature = float(given_temperature)

                # Stops the input prompt from repeating.
                successful = True
            # Catches error exception for when the given value is invalid.
            except ValueError:
                # Informs the user what values are valid.
                print(
                    "\nThe given temperature value must be a valid float or int ONLY."
                )

                # Repeats the input prompt.
                successful = False

        # Gets the thermodynamic values of the reactants and products.
        reactant_standard_heats, reactant_standard_entropies = (
            self.get_thermodynamic_values(self.reactants, "Reactant")
        )
        product_standard_heats, product_standard_entropies = (
            self.get_thermodynamic_values(self.products, "Product")
        )

        # Creates or opens a file to store where the results will be contained.
        file = open(file_path, "w+", encoding="utf-8")

        # Writes to the file the given chemical equation.
        file.write(f"The given chemical equation:\n{self.chemical_equation}")

        # Writes to the file the given temperature and the table of givens.
        file.write("\n\nGivens:")
        file.write(f"\nTemperature: {temperature} K")
        file.write("\n" + self.givens_table.get_string())  # type: ignore

        # Writes to the file the required values to be calculated.
        file.write("\n\nRequired:\nΔH°rxn, ΔS°rxn, and ΔS°univ")

        # Closes the created or opened file.
        file.close()

        # Calculates for the required values.
        standard_enthalpy = self.calculate_rxn_enthalpy(
            product_standard_heats, reactant_standard_heats
        )
        standard_entropy_change = self.calculate_rxn_entropy(
            product_standard_entropies, reactant_standard_entropies
        )
        surrounding_entropy = self.calculate_surr_entropy(
            standard_entropy_change,
            temperature,
        )
        universal_entropy = self.calculate_univ_entropy(
            standard_entropy_change, surrounding_entropy
        )

        # Opens the file where the results are contained.
        file = open(file_path, "a+", encoding="utf-8")

        # Writes to the file the final answers to the required values.
        file.write("\n\nThe FINAL ANSWERS:")
        file.write(f"\nΔH°rxn = {standard_enthalpy} J/K × mol")
        file.write(f"\nΔS°rxn = {standard_entropy_change} J/K × mol")
        file.write(f"\nΔS°univ = {universal_entropy} J/K × mol")

        # Closes the opened file.
        file.close()

        # Reopens the file where the results are contained for reading.
        file = open(file_path, "r", encoding="utf-8")

        # Reads the results contained in the file.
        file_contents = file.read()

        # Prints out the results.
        print("\n" + file_contents)

        # Closes the opened file.
        file.close()


class ChemicalEquilibrium:
    # Creates a new ASCII table to store the givens.
    givens_table = PrettyTable(
        ["Formula", "Type", "Number", "Molarity (M)", "Pressure (atm)"]
    )

    def __init__(self: Self, raw_chemical_equation: str) -> None:
        """Constructs a new instance of the ChemicalEquilibrium class,

        Args:
            self (Self): The class of the calculator.
            raw_chemical_equation (str): The chemical equation.
        """
        # Separates the chemical equation into its reactants and products.
        chemical_equation = raw_chemical_equation.replace(" ", "").split("<=>")
        # Formats and splits the reactants and products into individual formulas.
        self.reactants = list(map(subscript, chemical_equation[0].split("+")))
        self.products = list(map(subscript, chemical_equation[1].split("+")))
        # Rejoins the reactants and products into a stylized chemical equation.
        self.chemical_equation = " ⇌ ".join(
            [" + ".join(self.reactants), " + ".join(self.products)]
        )

    def get_values(
        self: Self, formulas: list[str], formula_type: str, type: str
    ) -> tuple[list[float], list[int]]:
        """Gets the molarity values of chemical formulas through inputs from the user.

        Args:
            self (Self): The class of the calculator.
            formulas (list[str]): The chemical formulas.
            formula_type (str): The type of the thermodynamic values.
            type (str): The type of equilibrium constant

        Returns:
            tuple[list[float], list[int]]: The user given values.
        """
        # Create lists to store the values and coefficients of the given formulas.
        values: list[float] = []
        coefficients: list[int] = []

        for formula in formulas:
            # Resets the variable to allow for looping if inputted value is invalid
            successful = False

            # Loops the input prompt until the user gives a valid value.
            while not successful:
                # Tries to prompt the user to input the values.
                try:
                    if formula[0].isdigit():
                        # Separates the parts of the given chemical formula into its coefficient and main.
                        coefficient = int(formula[0])
                        main_formula = formula[1:]

                        # Asks the user for the value and exponentiates it to the coefficient.
                        given_value = float(
                            input(
                                f'\nPlease input the {"MOLARITY (M)" if type == "c" else "PRESSURE (atm)"} of "{main_formula}" (float or int ONLY): '
                            )
                        )
                        value = given_value**coefficient

                        # Adds the given values to the table of givens
                        self.givens_table.add_row(  # type: ignore
                            [
                                main_formula,
                                formula_type,
                                coefficient,
                                f"{given_value} M" if type == "c" else "N/A",
                                f"{given_value} atm" if type == "p" else "N/A",
                            ]
                        )

                        coefficients.append(coefficient)
                    else:
                        # Asks the user for the value.
                        value = given_value = float(
                            input(
                                f'\nPlease input the {"MOLARITY (M)" if type == "c" else "PRESSURE (atm)"} of "{formula}" (float or int ONLY): '
                            )
                        )

                        # Adds the given values to the table of givens
                        self.givens_table.add_row(  # type: ignore
                            [
                                formula,
                                formula_type,
                                1,
                                f"{given_value} M" if type == "c" else "N/A",
                                f"{given_value} atm" if type == "p" else "N/A",
                            ]
                        )

                        coefficients.append(1)

                    # Adds the given values to the list that will be returned by the function
                    values.append(value)

                    # Stops the input prompt from repeating.
                    successful = True
                # Catches error exception for when the given value is invalid.
                except ValueError:
                    # Informs the user what values are valid.
                    print("\nThe given value must be a valid float or int ONLY.")

                    # Repeats the input prompt.
                    successful = False

        # Returns the list of given values and formula coefficients.
        return values, coefficients

    def calculate_m(
        self: Self,
        product_molarities: list[float],
        product_moles: list[int],
        reactant_molarities: list[float],
        reactant_moles: list[int],
    ) -> tuple[float, float]:
        """Calculates the equilibrium constant of the molarities of the reactants and products.

        Args:
            self (Self): The class of the calculator.
            product_molarities (list[float]): The molarities of the products.
            product_moles (list[int]): The moles of the products.
            reactant_molarities (list[float]): The molarities of the reactants.
            reactant_moles (list[int]):The moles of the products.

        Returns:
            float: The calculated equilibrium constant.
        """
        # Opens the file where the results are contained.
        file = open(file_path, "a", encoding="utf-8")

        # Writes to the file the equation that will be used to calculate the equilibrium constant.
        file.write("\n\nEquation for the equilibrium constant (Kc):\nKc = [B]ᵇ / [A]ᵃ")

        # Writes to the file the solution for calculating the equilibrium constant.
        file.write(
            f'\n\nSolution for the equilibrium constant (Kc) of "{self.chemical_equation}":'
        )
        file.write(
            f"\nKc = ({' × '.join(map(str, product_molarities))}) / ({' × '.join(map(str, reactant_molarities))})"
        )

        # Calculates the equilibrium constant.
        kc = reduce(lambda x, y: x * y, product_molarities) / reduce(
            lambda x, y: x * y, reactant_molarities
        )

        # Writes to the file the calculated equilibrium constant.
        file.write(f"\nKc = {kc}")

        file.write(f"\nKc = 1 / {kc}")

        # Calculates the backwards equilibrium constant.
        backward_kc = 1 / kc

        # Writes to the file the calculated backwards equilibrium constant.
        file.write(f"\nKc of backward = {backward_kc}")

        # Checks if the equilibrium constant is greater than zero then writes to the file that forward is favored.
        if kc > 0:
            file.write(
                "\n\nBecause the equilibrium constant (Kc) is greater than 0, FORWARD is more favored."
            )
        # Checks if the equilibrium constant is less than zero then writes to the file that backward is favored.
        elif kc < 0:
            file.write(
                "\n\nBecause the equilibrium constant (Kc) is less than 0, BACKWARD is more favored."
            )
        # Writes to the file that the reaction is at equilibrium.
        else:
            file.write(
                "\n\nBecause the equilibrium constant (Kc) is equal to 0, the reaction would not proceed and no products formed."
            )

        # Writes to the file the equation that will be used to calculate the equilibrium constant.
        file.write(
            "\n\nEquation for the equilibrium constant (Kp):\nKp = Kc(0.0821 × T)ⁿ"
        )

        # Writes to the file the solution for calculating the equilibrium constant.
        file.write(
            f'\n\nSolution for the equilibrium constant (Kp) of "{self.chemical_equation}":'
        )
        file.write("\nT = 800 C + 273.15\nT = 1073.15 K")
        file.write(f"\nΔn = {fsum(product_moles)} - {fsum(reactant_moles)}")

        # Calculates the number of moles.
        moles = fsum(product_moles) - fsum(reactant_moles)

        file.write(f"\nΔn = {moles} mole/s")
        file.write(f"\nKp = {kc}(0.0821 × 1073.15 K){superscript(str(int(moles)))}")

        # Calculates the equilibrium constant.
        kp = 0.0821 * 1073.15
        kp **= moles
        kp *= kc

        # Writes to the file the calculated equilibrium constant.
        file.write(f"\nKp = {kp}")

        # Closes the opened file.
        file.close()

        # Returns the calculated equilibrium constant.
        return kc, kp

    def calculate_p(
        self: Self,
        product_pressures: list[float],
        product_moles: list[int],
        reactant_pressures: list[float],
        reactant_moles: list[int],
    ) -> tuple[float, float]:
        """Calculates the equilibrium constant of the pressures of the reactants and products.

        Args:
            self (Self): The class of the calculator.
            product_pressures (list[float]): The pressures of the products.
            product_moles (list[int]): The moles of the products.
            reactant_pressures (list[float]): The pressures of the reactants.
            reactant_moles (list[int]):The moles of the products.

        Returns:
            float: The calculated equilibrium constant.
        """
        # Opens the file where the results are contained.
        file = open(file_path, "a", encoding="utf-8")

        # Writes to the file the equation that will be used to calculate the equilibrium constant.
        file.write(
            "\n\nEquation for the equilibrium constant (Kp):\nKp = [Pb]ᵇ / [Pa]ᵃ"
        )

        # Writes to the file the solution for calculating the equilibrium constant.
        file.write(
            f'\n\nSolution for the equilibrium constant (Kp) of "{self.chemical_equation}":'
        )
        file.write(
            f"\nKp = ({' × '.join(map(str, product_pressures))}) / ({' × '.join(map(str, reactant_pressures))})"
        )

        # Calculates the equilibrium constant.
        kp = reduce(lambda x, y: x * y, product_pressures) / reduce(
            lambda x, y: x * y, reactant_pressures
        )

        # Writes to the file the calculated equilibrium constant.
        file.write(f"\nKp = {kp}")

        # Writes to the file the equation that will be used to calculate the equilibrium constant.
        file.write(
            "\n\nEquation for the equilibrium constant (Kc):\nKc = Kp / (0.0821 × T)ⁿ"
        )

        # Writes to the file the solution for calculating the equilibrium constant.
        file.write(
            f'\n\nSolution for the equilibrium constant (Kc) of "{self.chemical_equation}":'
        )
        file.write("\nT = 800 C + 273.15\nT = 1073.15 K")
        file.write(f"\nΔn = {fsum(product_moles)} - {fsum(reactant_moles)}")

        # Calculates the number of moles.
        moles = fsum(product_moles) - fsum(reactant_moles)

        file.write(f"\nΔn = {moles} mole/s")
        file.write(f"\nKc = {kp} / (0.0821 × 1073.15 K){superscript(str(int(moles)))}")

        # Calculates the equilibrium constant.
        kc = 0.0821 * 1073.15
        kc **= moles
        kc = kp / kc

        # Writes to the file the calculated equilibrium constant.
        file.write(f"\nKc = {kc}")

        file.write(f"\nKc = 1 / {kc}")

        # Calculates the backwards equilibrium constant.
        backward_kc = 1 / kc

        # Writes to the file the calculated backwards equilibrium constant.
        file.write(f"\nKc of backward = {backward_kc}")

        # Checks if the equilibrium constant is greater than zero then writes to the file that forward is favored.
        if kc > 0:
            file.write(
                "\n\nBecause the equilibrium constant (Kc) is greater than 0, FORWARD is more favored."
            )
        # Checks if the equilibrium constant is less than zero then writes to the file that backward is favored.
        elif kc < 0:
            file.write(
                "\n\nBecause the equilibrium constant (Kc) is less than 0, BACKWARD is more favored."
            )
        # Writes to the file that the reaction is at equilibrium.
        else:
            file.write(
                "\n\nBecause the equilibrium constant (Kc) is equal to 0, the reaction would not proceed and no products formed."
            )

        # Closes the opened file.
        file.close()

        # Returns the calculated equilibrium constant.
        return kc, kp

    def is_valid(self: Self, formula: str) -> bool:
        """Checks whether a formula is valid or not.

        Args:
            self (Self): The class of the calculator.
            formula (str): The chemical formula.

        Returns:
            Boolean: Whether a formula is valid or not.
        """
        # Creates a list of chemical states to disregard.
        disregard = ["s)", "l)"]
        # Takes the formula's state.
        formula_state = formula.split("(")[1]

        # Checks if the formula should be disregarded.
        if formula_state in disregard:
            return False
        else:
            return True

    def calculate(self: Self) -> None:
        """Runs the chemical equilibrium calculator.

        Args:
            self (Self): The class of the calculator.
        """
        # Filters out the solid and liquid formulas.
        reactants = list(filter(self.is_valid, self.reactants))
        products = list(filter(self.is_valid, self.products))

        # Creates or opens a file to store where the results will be contained.
        file = open(file_path, "w+", encoding="utf-8")

        # Writes to the file the given chemical equation.
        file.write(f"The given chemical equation:\n{self.chemical_equation}")

        # Resets the variable to allow for looping if inputted value is invalid
        successful = False

        # Loops the input prompt until the user gives a valid "Y" or "n" value.
        while not successful:
            # Asks the user whether they would be using pressure or not.
            equilibrium_type = input(
                "Please input if you would be using pressure (Y/n): "
            )

            # Checks whether the given input is a "Y" value.
            if equilibrium_type == "Y":
                # Replaces the "Y" fot yes with the "p" for pressure.
                equilibrium_type = "p"

                # Stops the input prompt from repeating.
                successful = True
            # Checks whether the given input is a "n" value.
            elif equilibrium_type == "n":
                # Replaces the "n" fot no with the "c" for moles.
                equilibrium_type = "c"

                # Stops the input prompt from repeating.
                successful = True
            else:
                # Informs the user what values are valid.
                print('\nThe given value must be either "Y" or "n" only.')

                # Repeats the input prompt.
                successful = False
                continue

        # Gets the values of the reactants and products.
        reactant_values, reactant_moles = self.get_values(
            reactants,
            "Reactants",
            equilibrium_type,  # type: ignore
        )
        product_values, product_moles = self.get_values(
            products,
            "Products",
            equilibrium_type,  # type: ignore
        )

        # Writes to the file the table of givens.
        file.write("\n\nGivens:")
        file.write("\n" + self.givens_table.get_string())  # type: ignore

        # Checks if the values are pressures.
        if equilibrium_type == "p":  # type: ignore
            # Writes to the file the required values to be calculated.
            file.write("\n\nRequired:\nKp and Kc")

            # Closes the created or opened file.
            file.close()

            # Calculates for the required values.
            kc, kp = self.calculate_p(
                product_values,
                product_moles,
                reactant_values,
                reactant_moles,
            )

            # Opens the file where the results are contained.
            file = open(file_path, "a+", encoding="utf-8")

            # Writes to the file the final answers to the required values.
            file.write("\n\nThe FINAL ANSWERS:")
            file.write(f"\nKp = {kp}")
            file.write(f"\nKc = {kc}")
        else:
            # Writes to the file the required values to be calculated.
            file.write("\n\nRequired:\nKc and Kp")

            # Closes the created or opened file.
            file.close()

            # Calculates for the required values.
            kc, kp = self.calculate_m(
                product_values,
                product_moles,
                reactant_values,
                reactant_moles,
            )

            # Opens the file where the results are contained.
            file = open(file_path, "a+", encoding="utf-8")

            # Writes to the file the final answers to the required values.
            file.write("\n\nThe FINAL ANSWERS:")
            file.write(f"\nKc = {kc}")
            file.write(f"\nKp = {kp}")

        # Closes the opened file.
        file.close()

        # Reopens the file where the results are contained for reading.
        file = open(file_path, "r", encoding="utf-8")

        # Reads the results contained in the file.
        file_contents = file.read()

        # Prints out the results.
        print("\n" + file_contents)

        # Closes the opened file.
        file.close()


class AcidsAndBases:
    def calculate_concentration(self: Self, p: float, type: str) -> float:
        """Calculates the ion concentration of the potential.

        Args:
            self (Self): The class of the calculator.
            p (float): The potential.
            type (str): The type of ion concentration.

        Returns:
            float: The calculated ion concentration.
        """
        # Opens the file where the results are contained.
        file = open(file_path, "a", encoding="utf-8")

        if type == "H":
            # Writes to the file the equation that will be used to calculate the hydrogen ion concentration.
            file.write("\n\nEquation for the hydrogen ion (H⁺):\n[H⁺] = 10⁻ᵖᴴ")

            # Writes to the file the solution for calculating the hydrogen ion concentration.
            file.write(f'\n\nSolution for the  hydrogen ion (H⁺) of "pH = {p}":')
            file.write(f"\n[H⁺] = 10⁻{superscript(str(p))}")

            # Calculates the hydrogen ion concentration.
            h = p * -1
            h = 10**h

            # Writes to the file the calculated hydrogen ion concentration.
            file.write(f"\n[H⁺] = {h} M")
        else:
            # Writes to the file the equation that will be used to calculate the hydroxide ion concentration.
            file.write("\n\nEquation for the hydroxide ion (OH⁻):\n[OH⁻] = 10⁻ᵖᴼᴴ")

            # Writes to the file the solution for calculating the hydroxide ion concentration.
            file.write(f'\n\nSolution for the  hydroxide ion (OH⁻) of "pOH = {p}":')
            file.write(f"\n[OH⁻] = 10⁻{superscript(str(p))}")

            # Calculates the hydroxide ion concentration.
            h = p * -13
            h = 10**h

            # Writes to the file the calculated hydroxide ion concentration.
            file.write(f"\n[OH⁻] = {h} M")

        # Closes the opened file.
        file.close()

        # Returns the calculated ion concentration.
        return h

    def calculate_potential(self: Self, h: float, type: str) -> float:
        """Calculates the potential of the ion concentration.

        Args:
            self (Self): The class of the calculator.
            h (float): The ion concentration.
            type (str): The type of potential.

        Returns:
            float: The calculated potential.
        """
        # Opens the file where the results are contained.
        file = open(file_path, "a", encoding="utf-8")

        if type == "H":
            # Writes to the file the equation that will be used to calculate the potential of hydrogen ion.
            file.write("\n\nEquation for the hydrogen ion (H⁺):\npH = -log[H⁺]")

            # Writes to the file the solution for calculating the potential of hydrogen ion.
            file.write(f'\n\nSolution for the  hydrogen ion (H⁺) of "[H⁺] = {h}":')
            file.write(f"\npH = -log({h})")

            # Calculates the potential of hydrogen ion.
            p = log10(h) * -1

            # Writes to the file the calculated potential of hydrogen ion.
            file.write(f"\npH = {p}")
        else:
            # Writes to the file the equation that will be used to calculate the potential of hydroxide ion.
            file.write("\n\nEquation for the hydroxide ion (OH⁻):\npOH = -log[OH⁻]")

            # Writes to the file the solution for calculating the potential of hydroxide ion.
            file.write(f'\n\nSolution for the  hydroxide ion (OH⁻) of "[OH⁻] = {h}":')
            file.write(f"\npOH = -log({h})")

            # Calculates the potential of hydroxide ion.
            p = log10(h) * -1

            # Writes to the file the calculated potential of hydroxide ion.
            file.write(f"\npOH = {p}")

        # Closes the opened file.
        file.close()

        # Returns the calculated potential of ion concentration.
        return p

    def convert_potential(self: Self, p: float) -> float:
        """Converts the potential of hydrogen ion to the potential of hydroxide ion and vice-versa.

        Args:
            self (Self): The class of the calculator.
            p (float): The potential.

        Returns:
            float: The converted potential.
        """
        # Calculates the other potential and returns it.

        # Opens the file where the results are contained.
        file = open(file_path, "a", encoding="utf-8")

        # Calculates the potential conversion.
        p = 14 - p

        # Closes the opened file.
        file.close()

        # Returns the converted potential.
        return p

    def calculate(self: Self) -> None:
        """Runs the acids and bases calculator.

        Args:
            self (Self): The class of the calculator.
        """
        # Creates or opens a file to store where the results will be contained.
        file = open(file_path, "w+", encoding="utf-8")

        # Resets the variable to allow for looping if inputted value is invalid
        successful = False

        # Loops the input prompt until the user gives a valid value.
        while not successful:
            # Asks the user what value they have.
            required = input("\nPlease input what value you have (pH/H/pOH/OH): ")

            # Checks whether the given input is a "pH" value.
            if required == "pH":
                # Tries to prompt the user to input the molarity values.
                try:
                    ph = float(input("\nPlease input the pH (float or int ONLY): "))

                    # Writes to the file the given pH.
                    file.write("Givens:")
                    file.write(f"\npH = {ph}")

                    # Closes the created or opened file.
                    file.close()

                    h = self.calculate_concentration(ph, "H")
                    poh = self.convert_potential(ph)
                    oh = self.calculate_concentration(poh, "OH")

                    # Opens the file where the results are contained.
                    file = open(file_path, "a+", encoding="utf-8")

                    # Writes to the file the final answers to the required values.
                    file.write("\n\nThe FINAL ANSWERS:")
                    file.write(f"\n[H⁺] = {h}")
                    file.write(f"\npOH = {poh}")
                    file.write(f"\n[OH⁻] = {oh}")

                    # Stops the input prompt from repeating.
                    successful = True
                # Catches error exception for when the given value is invalid.
                except ValueError:
                    # Informs the user what values are valid.
                    print("\nThe given value must be a valid float or int ONLY.")

                    # Repeats the input prompt.
                    successful = False
            # Checks whether the given input is a "H" value.
            elif required == "H":
                # Tries to prompt the user to input the molarity values.
                try:
                    h = float(input("\nPlease input the H⁺ (float or int ONLY): "))

                    # Writes to the file the given pH.
                    file.write("Givens:")
                    file.write(f"\n[H⁺] = {h}")

                    # Closes the created or opened file.
                    file.close()

                    ph = self.calculate_potential(h, "H")
                    poh = self.convert_potential(ph)
                    oh = self.calculate_concentration(poh, "OH")

                    # Opens the file where the results are contained.
                    file = open(file_path, "a+", encoding="utf-8")

                    # Writes to the file the final answers to the required values.
                    file.write("\n\nThe FINAL ANSWERS:")
                    file.write(f"\npH = {ph}")
                    file.write(f"\npOH = {poh}")
                    file.write(f"\n[OH⁻] = {oh}")

                    # Stops the input prompt from repeating.
                    successful = True
                # Catches error exception for when the given value is invalid.
                except ValueError:
                    # Informs the user what values are valid.
                    print("\nThe given value must be a valid float or int ONLY.")

                    # Repeats the input prompt.
                    successful = False
            # Checks whether the given input is a "pOH" value.
            elif required == "pOH":
                # Tries to prompt the user to input the molarity values.
                try:
                    poh = float(input("\nPlease input the pOH (float or int ONLY): "))

                    # Writes to the file the given pH.
                    file.write("Givens:")
                    file.write(f"\npOH = {poh}")

                    # Closes the created or opened file.
                    file.close()

                    oh = self.calculate_concentration(poh, "OH")
                    ph = self.convert_potential(poh)
                    h = self.calculate_concentration(ph, "H")

                    # Opens the file where the results are contained.
                    file = open(file_path, "a+", encoding="utf-8")

                    # Writes to the file the final answers to the required values.
                    file.write("\n\nThe FINAL ANSWERS:")
                    file.write(f"\n[OH⁻] = {oh}")
                    file.write(f"\npH = {ph}")
                    file.write(f"\n[H⁺] = {h}")

                    # Stops the input prompt from repeating.
                    successful = True
                # Catches error exception for when the given value is invalid.
                except ValueError:
                    # Informs the user what values are valid.
                    print("\nThe given value must be a valid float or int ONLY.")

                    # Repeats the input prompt.
                    successful = False
            # Checks whether the given input is a ")H" value.
            elif required == "OH":
                # Tries to prompt the user to input the molarity values.
                try:
                    oh = float(input("\nPlease input the OH⁻ (float or int ONLY): "))

                    # Writes to the file the given pH.
                    file.write("Givens:")
                    file.write(f"\n[OH⁻] = {oh}")

                    # Closes the created or opened file.
                    file.close()

                    poh = self.calculate_potential(oh, "OH")
                    ph = self.convert_potential(poh)
                    h = self.calculate_concentration(ph, "H")

                    # Opens the file where the results are contained.
                    file = open(file_path, "a+", encoding="utf-8")

                    # Writes to the file the final answers to the required values.
                    file.write("\n\nThe FINAL ANSWERS:")
                    file.write(f"\npOH = {poh}")
                    file.write(f"\npH = {ph}")
                    file.write(f"\n[H⁺] = {h}")

                    # Stops the input prompt from repeating.
                    successful = True
                # Catches error exception for when the given value is invalid.
                except ValueError:
                    # Informs the user what values are valid.
                    print("\nThe given value must be a valid float or int ONLY.")

                    # Repeats the input prompt.
                    successful = False
            else:
                # Informs the user what values are valid.
                print('\nThe given value must be either "Y" or "n" only.')

                # Repeats the input prompt.
                successful = False
                continue

        # Reopens the file where the results are contained for reading.
        file = open(file_path, "r", encoding="utf-8")

        # Reads the results contained in the file.
        file_contents = file.read()

        # Prints out the results.
        print("\n" + file_contents)

        # Closes the opened file.
        file.close()


# Resets the variable to allow for looping if inputted value is invalid
successful = False

# Loops the input prompt until the user gives a valid chemical equation and temperature value.
while not successful:
    # Gives the user the instructions for providing the topic.
    print("""
Please input the corresponding integer for the General Chemistry Topic that you need this calculator for.
1.) Thermochemistry,
2.) Chemical Kinetics,
3.) Chemical Equilibrium, 
4.) Acids and Bases,
""")
    # Asks the user for the topic.
    topic = input("Topic: ")

    # Checks whether the given input is a "1" value.
    if topic == "1":
        # Resets the variable to allow for looping if inputted value is invalid
        successful = False

        # Loops the input prompt until the user gives a valid chemical equation.
        while not successful:
            # Gives the user the instructions for providing the chemical equation to use.
            print(
                "\nPlease input your chemical equation, this is case sensitive.\ne.g. N2 + 3H2 -> 2NH3"
            )
            # Asks the user for the chemical equation to use.
            given_chemical_equation = input("\nChemical Equation: ")

            # Checks if the there was no chemical equation given.
            if given_chemical_equation == "":
                # Informs the user that they need to provide a chemical equation.
                print("\nYou NEED provide a chemical equation.")

                # Repeats the input prompt.
                successful = False

            # Tries to run the chemical equilibrium calculator
            try:
                Thermochemistry(given_chemical_equation).calculate()  # type: ignore

                # Stops the input prompt from repeating.
                successful = True

            except IndexError:
                # Informs the user that their chemical equation does not follow the format.
                print("\nThe chemical equation MUST follow the necessary format.")

                # Repeats the input prompt.
                successful = False

        # Stops the input prompt from repeating.
        successful = True
    # Checks whether the given input is a "3" value.
    elif topic == "3":
        # Resets the variable to allow for looping if inputted value is invalid
        successful = False

        # Loops the input prompt until the user gives a valid chemical equation.
        while not successful:
            # Gives the user the instructions for providing the chemical equation to use.
            print(
                "\nPlease input your chemical equation, this is case sensitive.\ne.g. CaCO3(s) <=> CaO(s) + CO2(g)"
            )
            # Asks the user for the chemical equation to use.
            given_chemical_equation = input("\nChemical Equation: ")

            # Checks if the there was no chemical equation given.
            if given_chemical_equation == "":
                # Informs the user that they need to provide a chemical equation.
                print("\nYou NEED provide a chemical equation.")

                # Repeats the input prompt.
                successful = False

            # Tries to run the chemical equilibrium calculator
            try:
                ChemicalEquilibrium(given_chemical_equation).calculate()  # type: ignore

                # Stops the input prompt from repeating.
                successful = True
            except IndexError:
                # Informs the user that their chemical equation does not follow the format.
                print("\nThe chemical equation MUST follow the necessary format.")

                # Repeats the input prompt.
                successful = False

        # Stops the input prompt from repeating.
        successful = True
    # Checks whether the given input is a "4" value.
    elif topic == "4":
        # Runs the chemical equilibrium calculator
        AcidsAndBases().calculate()  # type: ignore

        # Stops the input prompt from repeating.
        successful = True
    else:
        # Informs the user what values are valid.
        print('\nThe given value must be "1", "2", "3", or "4" only.')

        # Repeats the input prompt.
        successful = False
        continue

# Resets the variable to allow for looping if inputted value is invalid
successful = False

# Loops the input prompt until the user gives a valid "Y" or "n" value.
while not successful:
    # Asks the user whether they want to keep the generated file or not.
    delete_file = input(
        "Do you want to keep the generated file containing the results? (Y/n): "
    )

    # Checks whether the given input is a "n" value.
    if delete_file == "n":
        # Deletes the generated file where the results are contained.
        remove(file_path)

        # Prints out a confirmation of the file's deletion.
        print("\nThe file was successfully deleted.")

        # Stops the input prompt from repeating.
        successful = True
    # Checks whether the given input is a "Y" value.
    elif delete_file == "Y":
        # Prints out the location of where the file was generated.
        print(f"\nThe file was generated at: {path.realpath(file_path)}\n")

        # Stops the input prompt from repeating.
        successful = True
    else:
        # Informs the user what values are valid.
        print('\nThe given value must be either "Y" or "n" only.')

        # Repeats the input prompt.
        successful = False
