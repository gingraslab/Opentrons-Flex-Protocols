#TODO: Serial dilution for standars

#standard slow down create more

metadata = {
    "protocolName": "NO QUICK TRANSFER AT ALL Single-plate BCA protocol with Quick Transfer",
    "author": "Nico To (modification of Sasha's original BCA protocol)",
    "description": "no quick transfer for protein addition step and standard creation. BCA for 1-24 samples.",
}
requirements = {"robotType": "Flex", "apiLevel": "2.20"}
import math
from opentrons import protocol_api

def get_vol_50ml_falcon(height):
    """
    Get's the volume of the liquid in the tube
    Height: height of liquid in tube in mm (start from tube bottom)
    Return: volume of liquid in tube in µl
    """
    volume = (1000 * (height - 9)) / 1.8
    return volume
def get_height_50ml_falcon(volume):
    """
    Get's the height of the liquid in the tube
    Volume: volume of liquid in tube in µl
    Return: hieght from bottom of tube in millimeters
    """
    height = (1.8 * (volume / 1000)) + 9
    return height
def get_height_15ml_falcon(volume):
    """
    Get's the height of the liquid in the tube
    Volume: volume of liquid in tube in ul
    Return: height in mm from the bottom of tube that pipette should go to
    """
    volume = volume / 1000
    if volume <= 1:  # cone part
        # print(-3.33*(volume**2)+15.45*volume+9.50)
        return -3.33 * (volume**2) + 15.45 * volume + 9.50 - 1  # −3.33x2+15.45x+9.50
    else:
        return 6.41667 * volume + 15.1667 - 5
def get_vol_15ml_falcon(height):
    """
    Get's the volume of the liquid in the tube
    Height: height in mm from the bottom of tube that pipette should go to
    Return: volume of liquid in tube in ul
    """
    if height <= 20.62:  # cone part
        volume = (((15.45 + math.sqrt(351.9225 - (13.32 * height)))) / 6.66) * 1000
        return volume
    else:
        volume = ((height - 10.1667) / 6.41667) * 1000
        return volume

def add_parameters(parameters):
    parameters.add_int(
        variable_name="number_samples",
        display_name="Number of Samples",
        description="Number of input samples.",
        default=10,
        minimum=1,
        maximum=40,
        unit="samples",
    )
    parameters.add_bool(
        variable_name="dulute_with_walt",
        display_name="Dilute using Walter",
        description="True: Dilute using Walter, False: Dilute manually",
        default=True,
    )
    parameters.add_int(
        variable_name="sample_vol",
        display_name="Amount of Sample",
        description="Amount of sample required in the duilution (only required if dilute_with_walt is True)",
        default=10,
        minimum=5,
        maximum=200,
        unit="ul",
    )
    parameters.add_int(
        variable_name="buffer_vol",
        display_name="Amount of Buffer",
        description="Amount of buffer required in the duilution (only required if dilute_with_walt is True)",
        default=90,
        minimum=5,
        maximum=200,
        unit="ul",
    )
    parameters.add_int(
        variable_name="working_sample_vol",
        display_name="Working Sample Volume",
        description="Volume of working sample (Volume for the samples in the flat well plate)",
        default=25,
        minimum=5,
        maximum=50,
        unit="ul",
    )
    parameters.add_int(
        variable_name="replication_mode",
        display_name="Replication Mode",
        choices=[
            {"display_name": "triplicate", "value": 3},
            {"display_name": "duplicate", "value": 2},
        ],
        default=3,
    )
    parameters.add_int(
        variable_name="incubation_time",
        display_name="Incubation Time",
        description="Incubation time",
        default=25,
        minimum=0,
        maximum=120,
        unit="minutes",
    )
    # parameters.add_int(
    #     variable_name="tip_type",
    #     display_name="types of tips",
    #     choices=[
    #         {"display_name": "1000ul", "value": 1000},
    #         {"display_name": "200ul", "value": 200},
    #     ],
    #     default=200,
    # )
    parameters.add_bool(
        variable_name="dry_run",
        display_name="Dry Run",
        description="Return tips and skip incubation (ignore this unless you are testing)",
        default=False,
    )


def run(protocol: protocol_api.ProtocolContext):
    replication_mode= protocol.params.replication_mode
    number_samples = protocol.params.number_samples
    is_dry_run = protocol.params.dry_run
    add_lid = True  # protocol.params.add_lid
    working_sample_vol = protocol.params.working_sample_vol
    pipette_max = 200-5

    # LOADING TIPS
    tips = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_200uL", slot)
        for slot in ["A3", "B3"]
    ]
    chute = protocol.load_waste_chute()

    def remove_tip(pipette):
        if is_dry_run:
            pipette.return_tip()
        else:
            pipette.drop_tip(chute)

    def find_aspirate_height(pip, source_well):
        lld_height = (
            pip.measure_liquid_height(source_well) - source_well.bottom().point.z
        )
        aspirate_height = max(lld_height - 5, 1)
        return aspirate_height

    # LOADING PIPETTES
    left_pipette = protocol.load_instrument(
        "flex_1channel_1000", "left", tip_racks=tips
    )
    right_pipette = protocol.load_instrument(
        "flex_8channel_1000", "right", tip_racks=tips
    )

    # LOADING LABWARE
    working_reagent_reservoir = protocol.load_labware("nest_12_reservoir_15ml", "D2")
    heatshaker = protocol.load_module("heaterShakerModuleV1", "D1")
    working_plate = protocol.load_labware("corning_96_wellplate_360ul_flat", "C2")
    reagent_stock = protocol.load_labware(
        "opentrons_10_tuberack_falcon_4x50ml_6x15ml_conical", "A1"
    )
    if add_lid:
        lid = protocol.load_labware(
            "opentrons_tough_pcr_auto_sealing_lid", location="C1"
        )
    sample_stock = protocol.load_labware(
        "opentrons_96_wellplate_200ul_pcr_full_skirt", "B2"
    )
    staging_slots = ["A4", "B4", "C4"]
    staging_racks = [
        protocol.load_labware("opentrons_flex_96_filtertiprack_1000uL", slot)
        for slot in staging_slots
    ]

    # REPLENISHING TIPS
    protein_addition_quick_transfer = False

    count = 0
    # DEFINING LIQUIDS
    bsa_stock = protocol.define_liquid(
        "BSA Stock", "BSA Stock from Pierce BCA Protein protocol ; 1.5mg/mL", "#FF6433"
    )
    bsa_rack = protocol.load_labware(
        "opentrons_24_tuberack_eppendorf_1.5ml_safelock_snapcap", "A2"
    )
    dye = protocol.define_liquid(
        "Reagent_A", "Reagent A for Working Reagent, will add 50 parts", "#FDF740"
    )
    # Reagent_B = protocol.define_liquid(
    #     "Reagent_B", "Reagent B for Working Reagent, will add 1 part", "#408AFD"
    # )
    water = protocol.define_liquid(
        "Diluent", "Diluent for standards, same as diluent in sample", "#D2E2FB"
    )
    sample = protocol.define_liquid(
        "sample",
        "Unknown Samples from CSV File",
        "#40FDF4",
    )

    empty_tube = protocol.define_liquid("empty", "Empty Tubes for Standards", "#D3D3D3")

    # LOADING LIQUIDS
    reagent_stock["A1"].load_liquid(water, 9000)
    bsa_rack["A1"].load_liquid(bsa_stock, 550)
    # reagent_stock["A3"].load_liquid(Reagent_A, 22000)
    # bsa_rack["D1"].load_liquid(Reagent_B, 1000)
    bsa_rack["B1"].load_liquid(empty_tube, 1)  # 1500 µg/mL
    bsa_rack["B2"].load_liquid(empty_tube, 1)  # 1000 µg/mL
    bsa_rack["B3"].load_liquid(empty_tube, 1)  # 750 µg/mL
    bsa_rack["B4"].load_liquid(empty_tube, 1)  # 500 µg/mL
    bsa_rack["B5"].load_liquid(empty_tube, 1)  # 250 µg/mL
    bsa_rack["B6"].load_liquid(empty_tube, 1)  # 125 µg/mL
    bsa_rack["C1"].load_liquid(empty_tube, 1)  # 25 µg/mL
    # bsa_rack["D6"].load_liquid(sample, 1)

    # dye_location = reagent_stock["A3"]
    dilutent_location = reagent_stock["A1"]
    sample_location = bsa_rack["D6"]
    bsa_stock_location = bsa_rack["A1"]
    
    heatshaker.open_labware_latch()
    #Diluting Sample
    diluted_sample_offset = 6
    

    if protocol.params.dulute_with_walt:
        left_pipette.pick_up_tip()
        vol_in_15_facon = get_vol_15ml_falcon(find_aspirate_height(left_pipette, dilutent_location))
        # num_transfers = math.ceil((number_samples*protocol.params.buffer_vol)/pipette_max)
        num_transfers = math.ceil((number_samples*protocol.params.buffer_vol)/(protocol.params.buffer_vol*math.floor(pipette_max/protocol.params.buffer_vol)))
        well_counter = 0
        col_num = replication_mode+1     # col num for the working_plate
        # print(num_transfers)
        for i in range (0, num_transfers):
            if i != num_transfers-1:    # not on last iteration
                aspirate_vol = pipette_max - pipette_max%protocol.params.buffer_vol
            else:
                aspirate_vol = (number_samples*protocol.params.buffer_vol)-(pipette_max - pipette_max%protocol.params.buffer_vol)*(num_transfers-1)
            if left_pipette.has_tip == False:
                left_pipette.pick_up_tip()
            left_pipette.blow_out(dilutent_location.top())
            left_pipette.aspirate(aspirate_vol+5, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 1)
            for x in range (0, math.floor(aspirate_vol/protocol.params.buffer_vol)):
                left_pipette.dispense(protocol.params.buffer_vol, sample_stock.wells()[well_counter + 48], 0.75)
                well_counter += 1
            # remove_tip(left_pipette)
            vol_in_15_facon-=aspirate_vol+5
        remove_tip(left_pipette)
        for i in range (0, math.ceil(number_samples/8)):
            right_pipette.pick_up_tip()
            right_pipette.aspirate(protocol.params.sample_vol, sample_stock['A' + str(i+1)].bottom(0.1), 0.5)
            right_pipette.dispense(protocol.params.sample_vol, sample_stock['A' + str(i+1+diluted_sample_offset)], 0.5)
            right_pipette.mix(3, protocol.params.sample_vol + protocol.params.buffer_vol-10, sample_stock['A' + str(i+1+diluted_sample_offset)], 0.5)
            # right_pipette.blow_out(sample_stock['A' + str(i+1+diluted_sample_offset)].top())
            right_pipette.touch_tip(sample_stock['A' + str(i+1+diluted_sample_offset)])
            if protein_addition_quick_transfer:
                right_pipette.aspirate(working_sample_vol*replication_mode+10, sample_stock['A' + str(i+1+diluted_sample_offset)],0.3)
                for x in range (0,replication_mode):
                    right_pipette.dispense(working_sample_vol, working_plate['A' + str(col_num)].bottom(0.5), 0.5)
                    # right_pipette.blow_out(working_plate['A' + str(col_num)].top())
                    col_num+=1
                remove_tip(right_pipette)
            else:
                for x in range (0,replication_mode):
                    if right_pipette.has_tip == False:
                        right_pipette.pick_up_tip()
                    right_pipette.aspirate(working_sample_vol, sample_stock['A' + str(i+1+diluted_sample_offset)],0.3)
                    right_pipette.dispense(working_sample_vol, working_plate['A' + str(col_num)].bottom(0.5), 0.5)
                    right_pipette.blow_out(working_plate['A' + str(col_num)].top())
                    col_num+=1
                    remove_tip(right_pipette)
            
            
    else:
        col_num = replication_mode+1
        if protein_addition_quick_transfer:
            for i in range (0, math.ceil(number_samples/8)):
                right_pipette.pick_up_tip()
                right_pipette.aspirate(working_sample_vol*replication_mode+10, sample_stock['A' + str(i+1)],0.5)
                for x in range (0,replication_mode):
                    right_pipette.dispense(working_sample_vol, working_plate['A' + str(col_num)].bottom(0.5), 0.5)
                    # right_pipette.blow_out(working_plate['A' + str(col_num)].top())
                    col_num+=1
                remove_tip(right_pipette)
        else:
            for i in range (0, math.ceil(number_samples/8)):
                # right_pipette.pick_up_tip()
                for x in range (0,replication_mode):
                    if right_pipette.has_tip == False:
                        right_pipette.pick_up_tip()
                    right_pipette.aspirate(working_sample_vol, sample_stock['A' + str(i+1)],0.5)
                    # right_pipette.aspirate(working_sample_vol, sample_stock['A' + str(i+1+diluted_sample_offset)],0.3)
                    right_pipette.dispense(working_sample_vol, working_plate['A' + str(col_num)].bottom(0.5), 0.5)
                    right_pipette.blow_out(working_plate['A' + str(col_num)].top())
                    col_num+=1
                    remove_tip(right_pipette)
                # remove_tip(right_pipette)

            
    def standard_loading(old, new):
        """
        old: well from sample stock
        new: row letter from sample plate
        """
        # left_pipette.pick_up_tip()
        if protein_addition_quick_transfer:
            left_pipette.aspirate(working_sample_vol*replication_mode+5, bsa_rack[old].bottom(1.5), 0.25)

            for i in range(1, replication_mode+1):  # A1,A2,A3
                left_pipette.dispense(working_sample_vol, working_plate[new + str(i)].bottom(0.1), 0.25)
        else:

            for i in range(1, replication_mode+1):  # A1,A2,A3
                if left_pipette.has_tip == False:
                    left_pipette.pick_up_tip()
                left_pipette.aspirate(working_sample_vol, bsa_rack[old].bottom(1.5), 0.25)
                left_pipette.dispense(working_sample_vol, working_plate[new + str(i)].bottom(0.1), 0.25)
                remove_tip(left_pipette)

    # Standard Preparation  FINISH LATER
    # standard_vol_per_tube = 500#working_sample_vol*replication_mode+50
    standard_vol_per_tube = working_sample_vol*replication_mode+60
    # dilutent_percentages = [0.25, 0.5, 0.625, 0.75, 0.875, 0.9375, 0.9875]
    dilutent_percentages = [1, 2/3, 1/2, 1/3, 1/6, 1/12, 1/60]
    buffer_vols =[]
    bsa_vols = []
    dilutent_pipette_vols = []
    serial_dilution_stock_tube = 3      # tubes 0 to 7
    amt_extra_in_1_over_12_tube=0
    for i in range(0, len(dilutent_percentages)):
        if standard_vol_per_tube*dilutent_percentages[i] < 5:
            #add extra to the serial_dilution_stock_tube duluted tube to make up for the 5ul minimum
            amt_extra_in_1_over_12_tube = (1/dilutent_percentages[serial_dilution_stock_tube])*standard_vol_per_tube*dilutent_percentages[i]
            amt_in_1_over_12_tube = standard_vol_per_tube + amt_extra_in_1_over_12_tube
            bsa_vols[serial_dilution_stock_tube] = amt_in_1_over_12_tube*dilutent_percentages[serial_dilution_stock_tube]
            buffer_vols[serial_dilution_stock_tube] = amt_in_1_over_12_tube*(1-dilutent_percentages[serial_dilution_stock_tube])
            bsa_vols.append(0)
            buffer_vols.append(standard_vol_per_tube-amt_extra_in_1_over_12_tube)
            # buffer_vols.append(standard_vol_per_tube*(1-dilutent_percentages[i]) - amt_extra_in_1_over_12_tube*(1-dilutent_percentages[i]))
            # print(standard_vol_per_tube*(1-dilutent_percentages[i]) - amt_extra_in_1_over_12_tube*(1-dilutent_percentages[i]))

        else:
            bsa_vols.append(standard_vol_per_tube*dilutent_percentages[i])
            buffer_vols.append(standard_vol_per_tube*(1-dilutent_percentages[i]))
    
    
    print(buffer_vols)
    print(bsa_vols)
    total_dilutent = 0
    for i in range(0, len(buffer_vols)):
        if total_dilutent + buffer_vols[i] < (pipette_max - 10):
            total_dilutent += buffer_vols[i]
            # print(total_dilutent)
            if i == len(buffer_vols) - 1:
                dilutent_pipette_vols.append(total_dilutent)
        else:
            dilutent_pipette_vols.append(total_dilutent)
            total_dilutent = buffer_vols[i]
            if i == len(buffer_vols) - 1:
                dilutent_pipette_vols.append(total_dilutent)
    # dilutent_pipette_vols.append(total_dilutent)
    
    # left_pipette.pick_up_tip()
    tube_spots = ["B1", "B2", "B3", "B4", "B5", "B6", "C1"]
    well_num = 0
    left_pipette.pick_up_tip()
    try:
        vol_in_15_facon
    except NameError:
        vol_in_15_facon = get_vol_15ml_falcon(find_aspirate_height(left_pipette, dilutent_location))
    if protein_addition_quick_transfer:
        for i in range(0, len(dilutent_pipette_vols)):
            vol_in_15_facon -= dilutent_pipette_vols[i]+10
            left_pipette.aspirate(
                dilutent_pipette_vols[i] + 10, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 0.5
            )
            amt_in_tip = dilutent_pipette_vols[i] + 10
            while amt_in_tip > 10:
                try:        # someties amt_in_tip is 10.0000000001
                    left_pipette.dispense(
                        buffer_vols[well_num],
                        bsa_rack[tube_spots[well_num]],
                        0.5,
                    )
                    amt_in_tip -= buffer_vols[well_num]
                except:
                    # print(amt_in_tip)
                    break

                well_num += 1
        remove_tip(left_pipette)
    else:
        rack_order = ["B1", "B2", "B3", "B4", "B5", "B6", "C1"]
        well_order = ["A", "B", "C", "D", "E", "F", "G"]
        for i in range (0, len(buffer_vols)):
            if left_pipette.has_tip == False:
                left_pipette.pick_up_tip()
            vol_in_15_facon -= buffer_vols[i]
            left_pipette.aspirate(
                buffer_vols[i], dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 0.5
            )
            left_pipette.dispense(buffer_vols[i],bsa_rack[tube_spots[i]])
            left_pipette.blow_out(bsa_rack[rack_order[i]].top())
            remove_tip(left_pipette)

    rack_order = ["B1", "B2", "B3", "B4", "B5", "B6", "C1"]
    well_order = ["A", "B", "C", "D", "E", "F", "G"]
    
    for i in range(0, len(bsa_vols)):
        left_pipette.pick_up_tip()
        if bsa_vols[i] == 0:
            left_pipette.aspirate(amt_extra_in_1_over_12_tube, bsa_rack[tube_spots[serial_dilution_stock_tube]], 0.5)
            print(amt_extra_in_1_over_12_tube)
        else:
            left_pipette.aspirate(
                bsa_vols[i],
                bsa_stock_location,
                0.5,
            )
        left_pipette.dispense(
            bsa_vols[i],
            bsa_rack[tube_spots[i]],
            0.5,
        )
        left_pipette.mix(4, standard_vol_per_tube - 10, bsa_rack[tube_spots[i]], 0.5)
        
        standard_loading(rack_order[i], well_order[i])
        # left_pipette.blow_out(bsa_rack[tube_spots[i]])
        if left_pipette.has_tip:
            remove_tip(left_pipette)

    # Vial H: Blank
    if protein_addition_quick_transfer:
        left_pipette.pick_up_tip()
        vol_in_15_facon-=working_sample_vol*replication_mode+5
        left_pipette.aspirate(working_sample_vol*replication_mode+5, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 0.25)
        for i in range(1, replication_mode+1):  # A1,A2,A3
            left_pipette.dispense(working_sample_vol, working_plate["H" + str(i)].bottom(0.1), 0.25)
        remove_tip(left_pipette)
    else:
        for i in range(1, replication_mode+1):  # A1,A2,A3
            left_pipette.pick_up_tip()
            vol_in_15_facon-=working_sample_vol
            left_pipette.aspirate(working_sample_vol, dilutent_location.bottom(get_height_15ml_falcon(vol_in_15_facon)), 0.25)
            left_pipette.dispense(working_sample_vol, working_plate["H" + str(i)].bottom(0.1), 0.25)
            remove_tip(left_pipette)
        
    # Adding Working Reagent to Plate
    num_columns = math.ceil(((math.ceil(number_samples / 8) * 8) * replication_mode) / 8) + replication_mode
    working_reagent_volume = 200
    right_pipette.pick_up_tip()
    working_reagent_volume_amt = num_columns*200*8 +1000#(working_reagent_volume*8*(math.ceil(number_samples/8)) + 1000)/8
    # Loading liquid for protocol setup
    for i in range (0, math.ceil(working_reagent_volume_amt/(10.5*1000))):
        if i ==math.ceil(working_reagent_volume_amt/(10.5*1000))-1:
            working_reagent_reservoir.wells()[i].load_liquid(dye, working_reagent_volume_amt-10500*i + 1500)
        else:
            working_reagent_reservoir.wells()[i].load_liquid(dye, 12000)

    for i in range (0, num_columns):
        working_reagent_volume_amt = working_reagent_volume_amt-(working_reagent_volume*8)
        # protocol.comment(str(working_reagent_volume_amt))
        # protocol.comment(str(math.ceil(working_reagent_volume_amt/(10.5*1000))))
        right_pipette.aspirate(working_reagent_volume, working_reagent_reservoir['A'+str(math.ceil(working_reagent_volume_amt/(10.5*1000)))],0.5)
        right_pipette.dispense(working_reagent_volume, working_plate["A"+str(i+1)].top(-1), rate=0.3)
        right_pipette.blow_out(working_plate["A"+str(i+1)].top(-1))
        right_pipette.blow_out(working_plate["A"+str(i+1)].top(-1))
    remove_tip(right_pipette)

    # Prep HeaterShaker
    heatshaker.open_labware_latch()
    # move labware with lid onto hs_mod
    if add_lid:
        protocol.deck.__delitem__("C2")
        new_working_plate = protocol.load_labware(
            "opentrons_96_wellplate_200ul_pcr_full_skirt", "C2"
        )
        protocol.move_labware(
            new_working_plate, new_location=heatshaker, use_gripper=True
        )
        heatshaker.close_labware_latch()
        heatshaker.open_labware_latch()
        new_working_plate.set_offset(x=0.00, y=0.00, z=30)
        protocol.move_labware(
            labware=lid, new_location=new_working_plate, use_gripper=True
        )
    else:
        protocol.move_labware(
            labware=working_plate, new_location=heatshaker, use_gripper=True
        )
    # protocol.pause("Place lid on well plate")
    heatshaker.close_labware_latch()
    heatshaker.set_and_wait_for_temperature(37)
    heatshaker.set_and_wait_for_shake_speed(400)

    # Shake For 30 Seconds
    if is_dry_run:
        protocol.delay(seconds=10)
    else:
        protocol.delay(minutes=0.5)
    heatshaker.deactivate_shaker()

    protocol.comment("\n---------------25 Minute Incubation----------------\n\n")
    if is_dry_run:
        protocol.delay(seconds=10)
    else:
        protocol.delay(minutes=protocol.params.incubation_time)  # SEND EMAIL AT 10 MINUTES

    # Deactivating Heatshaker
    heatshaker.deactivate_heater()
    heatshaker.open_labware_latch()
    if add_lid:
        protocol.move_labware(lid, "C1", use_gripper=True)
        protocol.move_labware(new_working_plate, "C2", use_gripper=True)
        heatshaker.close_labware_latch()
        # left_pipette.pick_up_tip()
        # left_pipette.pick_up_tip()
    else:
        protocol.move_labware(working_plate, "C2", use_gripper=True)
        heatshaker.close_labware_latch()
        
        # heatshaker.close_labware_latch()
