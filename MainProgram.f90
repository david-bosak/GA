! .#####...#####...........##...##...####...######...####...##..##.
! .##..##..##..##..........###.###..##..##....##....##..##..##..##.
! .##..##..#####...........##.#.##..######....##....##......######.
! .##..##..##..##....##....##...##..##..##....##....##..##..##..##.
! .#####...#####.....##....##...##..##..##....##.....####...##..##.
! .................................................................
! 
! Description:
!   Performance class containing all functions required to evaluate off-design performance of a gas-turbine system 
!   by embedding a genetic algorithm optimiser within it.
!   Should a decision be made to use a different method of optimisation, then this is the only class that needs to be modified.
!   All other dependent classes should not be affected.
!
! The process:
!     Opertiaon 1. Process inputs
!     Operation 2. Initialise files and folder structure
!     Operation 3. Generate initial population, update log file
!     Operation 4. Stochastic selection process
!     Operation 5. Cross over process
!     Operation 6. Mutation process
!     Operation 7. Check constraints
!     Operation 8. Write output, update log file
! 
! 	Version:	3.2.0
! 	Date:		24/03/2017
! 	Author:		David Bosak
! 	Changes:
! 				1. 	Fixed bugs in convergence criteria. some instances a brick can return a NaN local fitness function value, which
!                   correspond to insufficient temperature profile in Log Mean Temperature Difference calculations.
!                   It was found that instead of setting a constraint limit on the temperatures within
!                   individual bricks, it is more efficient to return and accept a NaN local fitness function and
!                   then address the issue on a global scale
!                
! 	Version:	3.0.0
! 	Date:		29/10/2016
! 	Author:		David Bosak
! 	Changes:
!               1. Improved the convergence. The system of non-linear equations can be reduced by taking away
!                   the Stodola like Law of Cones, which significantly improves the efficiency of the algorithm, reduces
!                   system complexity, and stabilises the progression of the solution across generations
!               2. Improved log functionality. The results are saved in three distinct logs - fitness.txt, genotypes.txt, phenotypes.txt
! 
! 	Version:	2.0.0
! 	Date:		16/12/2015
! 	Author:		David Bosak
! 	Changes:
! 				1.	Embedded genetic algorithm process within gas turbine off-design performance model
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program MainProgram
    
    
    USE DFLIB !library to create a directory or delete files
    USE VarDeclar
    
    implicit none
    CALL RANDOM_SEED !library to generate random number
    CALL ReadVar
    CALL ArrayInitialisation
    CALL ArrayAllocation
    CALL ZScoreTableInitialisation
    
    !initialise the Zscores array. 
    DummyReal=-4.09
    DO I=1, 819
        ZscoresArray(I,1) = DummyReal
        DummyReal = DummyReal + 0.01
    END DO
        

!-------------------------------------------------------- Process inputs, and initialise files and folders structure
!delete files and folders
    DO I=1, NoOfGenerations
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I))//"\Phenotypes.txt")
        !IF (result) THEN
        !    WRITE (*,*) I, ' Phenotypes.txt successfully deleted'
        !ELSE
        !    WRITE (*,*) I, ' Phenotypes.txt failed to delete'
        !END IF
        
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I))//"\Genotypes.txt")
        !IF (result) THEN
        !    WRITE (*,*) I, ' Phenotypes.txt successfully deleted'
        !ELSE
        !    WRITE (*,*) I, ' Phenotypes.txt failed to delete'
        !END IF
        
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I))//"\Fitness.txt")
        !IF (result) THEN
        !    WRITE (*,*) I, ' Phenotypes.txt successfully deleted'
        !ELSE
        !    WRITE (*,*) I, ' Phenotypes.txt failed to delete'
        !END IF
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt")
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_errors.txt")
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_statistics.txt")
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_full_genotypes.txt")
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_full_phenotypes.txt")
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_full_fitness.txt")
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_min_genotypes.txt")
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_min_phenotypes.txt")
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_min_fitness.txt")
        result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\NR_Log_Variables.txt")
        
        
        result = DelDirQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I)))
        
        
        !IF (result) THEN
        !    WRITE (*,*) I, ' folder successfully deleted'
        !ELSE
        !    WRITE (*,*) I, ' folder failed to delete'
        !END IF
        
    END DO
    !result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_NR.txt")
!Create folders
    DO I=1, NoOfGenerations
        result = MakeDirQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I)))
        !IF (result) THEN
        !    WRITE (*,*) I, ' folder successfully created'
        !ELSE
        !    WRITE (*,*) I, ' folder failed to create'
        !END IF
    END DO  
    
    !!!!!!!!!!! Delete all files in \variable folder
    NoOfBricksSteamPath_results_OD = -1
    NoOfBricksGasPath_results_OD = 0
    OPEN (150,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_results.txt",ACTION="READ")
        READ(150,*) DummyCharacter
        DummyCharacter = ""
        DO WHILE(DummyCharacter /= "Brick")
            READ(150,*) DummyCharacter
            NoOfBricksSteamPath_results_OD = NoOfBricksSteamPath_results_OD + 1
        END DO
        
        DO WHILE(.TRUE.)
            READ(150,*,IOSTAT=DummyInteger) DummyCharacter
            IF (DummyInteger.NE.0) THEN
                EXIT
            ELSE
                NoOfBricksGasPath_results_OD = NoOfBricksGasPath_results_OD + 1
            END IF
        END DO
    CLOSE(150)

    OPEN (150,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_results.txt",ACTION="READ")
        READ(150,*) DummyCharacter
        DO I=1, NoOfBricksSteamPath_results_OD
            READ(150,*)                         DummyCharacter, &
                                                DummyInteger
            result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(DummyInteger))//"_results_OD.txt")
            result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(DummyInteger))//"_results_OD_ff.txt")
            !!! Create a file to store fitness function values
            !OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(DummyInteger))//"_results_OD_ff.txt",ACTION="WRITE")
            !    WRITE(200,"(9999(G15.2,:))") "FF", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9"
            !CLOSE(200) 
        END DO
        
        !for gas path
        READ(150,*) DummyCharacter
        DO I=1, NoOfBricksGasPath_results_OD
            READ(150,*)                         DummyCharacter, &
                                                DummyInteger
            result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(DummyInteger))//"_results_OD.txt")
            result = DelFilesQQ(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(DummyInteger))//"_results_OD_ff.txt")
        END DO
    CLOSE(150)
    
    
                        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------- BEGIN OPTIMISATION ----------------------------------

!----------------------------------------------------- Generate Initial Population

        Print*, "Generating initial population" 
DO P=1, PopulationSize
    V = 0 !index used to check if fitness is NaN, if it is NaN then generate a new Individual.
    DO WHILE (V == 0)
            Print*, "Individual: ", P
                            OPEN (UNIT=200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_input.txt", ACTION="READ")  
                                !read stations    
                                    DO WHILE(DummyCharacter /= "start.")
                                        READ(200,*) DummyCharacter
                                    END DO !END, start. is found
                                    J = 0
                                    K = 0
                                    M = 0
                                    DO I=1,NoOfBricks_InputFile
                                        READ(200,*)     BrickName_InputFile(I), & !store names of bricks
                                                        DummyCharacter, &
                                                        StationNumberSteam_InputFile(I,1), & !store frst inlet steam station
                                                        StationNumberSteam_InputFile(I,2), & !store second inlet steam station
                                                        StationNumberSteam_InputFile(I,3), & !store first outlet steam station
                                                        StationNumberSteam_InputFile(I,4), & !store second outlet steam station
                                                        DummyCharacter, &
                                                        StationNumberGas_InputFile(I,1), & !store inlet gas station
                                                        StationNumberGas_InputFile(I,2)  !store outlet gas station
                
                                        IF (BrickName_InputFile(I) == "Econ") THEN
                                            Brick_Phenotypes(J+1) = BrickName_InputFile(I) !location reserved for Tg_out
                                            Brick_Phenotypes(J+2) = BrickName_InputFile(I) !location reserved for Ts_out
                                            Brick_Phenotypes(J+3) = BrickName_InputFile(I) !location reserved for Psout
                                            Brick_Phenotypes(J+4) = BrickName_InputFile(I) !location reserved for MFsout
                                            Brick_Phenotypes(J+5) = BrickName_InputFile(I) !location reserved for MFgout
                                            StationNumber_Phenotypes(J+1) = StationNumberGas_InputFile(I,2) !location reserved for Tg_out
                                            StationNumber_Phenotypes(J+2) = StationNumberSteam_InputFile(I,3) !location reserved for Ts_out
                                            StationNumber_Phenotypes(J+3) = StationNumberSteam_InputFile(I,3) !location reserved for Psout
                                            StationNumber_Phenotypes(J+4) = StationNumberSteam_InputFile(I,3) !location reserved for MFsout
                                            StationNumber_Phenotypes(J+5) = StationNumberGas_InputFile(I,2) !location reserved for MFgout
                                            !Variable index - 1=T, 2=P, 3=h, 4=S, 5=MS (temperature, pressure, enthalpy, entropy, flow)
                                            VariableType_Phenotypes(J+1) = 1 !array index reserved for Tg_out
                                            VariableType_Phenotypes(J+2) = 1 !array index reserved for Ts_out
                                            VariableType_Phenotypes(J+3) = 2 !array index reserved for Psout
                                            VariableType_Phenotypes(J+4) = 5 !array index reserved for MFsout
                                            VariableType_Phenotypes(J+5) = 5 !array index reserved for MFgout
                                            
                                            !Initialising in the first generation
                                            Phenotypes(J+4) = rand(10.0,100.0) !array index reserved for MFsout
                                            Phenotypes(J+5) = MFg_exhaustGT_OD !array index reserved for MFgout
                                            
                                            J = J + 5 !there are 5 variables in Econ
                    
                                        ELSE IF (BrickName_InputFile(I) == "Evap") THEN
                                            Brick_Phenotypes(J+1) = BrickName_InputFile(I) !location reserved for Tg_out
                                            Brick_Phenotypes(J+2) = BrickName_InputFile(I) !location reserved for Ts_out
                                            Brick_Phenotypes(J+3) = BrickName_InputFile(I) !location reserved for Psout
                                            Brick_Phenotypes(J+4) = BrickName_InputFile(I) !location reserved for MFsout
                                            Brick_Phenotypes(J+5) = BrickName_InputFile(I) !location reserved for MFgout
                                            Brick_Phenotypes(J+6) = BrickName_InputFile(I) !location reserved for Tsin
                                            Brick_Phenotypes(J+7) = BrickName_InputFile(I) !location reserved for Psin
                                            StationNumber_Phenotypes(J+1) = StationNumberGas_InputFile(I,2) !location reserved for Tg_out
                                            StationNumber_Phenotypes(J+2) = StationNumberSteam_InputFile(I,3) !location reserved for Ts_out
                                            StationNumber_Phenotypes(J+3) = StationNumberSteam_InputFile(I,3) !location reserved for Psout
                                            StationNumber_Phenotypes(J+4) = StationNumberSteam_InputFile(I,3) !location reserved for MFsout
                                            StationNumber_Phenotypes(J+5) = StationNumberGas_InputFile(I,2) !location reserved for MFgout
                                            StationNumber_Phenotypes(J+6) = StationNumberSteam_InputFile(I,1) !location reserved for Tsin
                                            StationNumber_Phenotypes(J+7) = StationNumberSteam_InputFile(I,1) !location reserved for Psin
                                            !Variable index - 1=T, 2=P, 3=h, 4=S, 5=MS (temperature, pressure, enthalpy, entropy, flow)
                                            VariableType_Phenotypes(J+1) = 1 !array index reserved for Tg_out
                                            VariableType_Phenotypes(J+2) = 1 !array index reserved for Ts_out
                                            VariableType_Phenotypes(J+3) = 2 !array index reserved for Psout
                                            VariableType_Phenotypes(J+4) = 5 !array index reserved for MFsout
                                            VariableType_Phenotypes(J+5) = 5 !array index reserved for MFgout
                                            VariableType_Phenotypes(J+6) = 1 !array index reserved for Tsin
                                            VariableType_Phenotypes(J+7) = 2 !array index reserved for Psin
                                            
                                            !Initialising in the first generation
                                            Phenotypes(J+4) = rand(10.0,100.0) !array index reserved for MFsout
                                            Phenotypes(J+5) = MFg_exhaustGT_OD !array index reserved for MFgout

                                            J = J + 7 !there are 7 variables in Evap
                    
                                        ELSE IF (BrickName_InputFile(I) == "SuHe") THEN
                                            Brick_Phenotypes(J+1) = BrickName_InputFile(I) !location reserved for Tg_out
                                            Brick_Phenotypes(J+2) = BrickName_InputFile(I) !location reserved for Ts_out
                                            Brick_Phenotypes(J+3) = BrickName_InputFile(I) !location reserved for Psout
                                            Brick_Phenotypes(J+4) = BrickName_InputFile(I) !location reserved for MFsout
                                            Brick_Phenotypes(J+5) = BrickName_InputFile(I) !location reserved for MFgout
                                            StationNumber_Phenotypes(J+1) = StationNumberGas_InputFile(I,2) !location reserved for Tg_out
                                            StationNumber_Phenotypes(J+2) = StationNumberSteam_InputFile(I,3) !location reserved for Ts_out
                                            StationNumber_Phenotypes(J+3) = StationNumberSteam_InputFile(I,3) !location reserved for Psout
                                            StationNumber_Phenotypes(J+4) = StationNumberSteam_InputFile(I,3) !location reserved for MFsout
                                            StationNumber_Phenotypes(J+5) = StationNumberGas_InputFile(I,2) !location reserved for MFgout
                                            !Variable index - 1=T, 2=P, 3=h, 4=S, 5=MS (temperature, pressure, enthalpy, entropy, flow)
                                            VariableType_Phenotypes(J+1) = 1 !array index reserved for Tg_out
                                            VariableType_Phenotypes(J+2) = 1 !array index reserved for Ts_out
                                            VariableType_Phenotypes(J+3) = 2 !array index reserved for Psout
                                            VariableType_Phenotypes(J+4) = 5 !array index reserved for MFsout
                                            VariableType_Phenotypes(J+5) = 5 !array index reserved for MFgout
                                            
                                            !Initialising in the first generation
                                            Phenotypes(J+4) = rand(10.0,100.0) !array index reserved for MFsout
                                            Phenotypes(J+5) = MFg_exhaustGT_OD !array index reserved for MFgout
                                            
                                            J = J + 5 !there are 5 variables in Suhe
                    
                                        ELSE IF (BrickName_InputFile(I) == "StTu") THEN
                                            Brick_Phenotypes(J+1) = BrickName_InputFile(I) !location reserved for MFsin
                                            StationNumber_Phenotypes(J+1) = StationNumberSteam_InputFile(I,1) !array index reserved for MFsin
                                            !Variable index - 1=T, 2=P, 3=h, 4=S, 5=MS (temperature, pressure, enthalpy, entropy, flow)
                                            VariableType_Phenotypes(J+1) = 5 !array index reserved for MFsin
                    
                                            !Initialising in the first generation
                                            Phenotypes(J+1) = rand(10.0,100.0) !array index reserved for MFsin

                                            J = J + 1 !there is 1 variables in StTu
                    
                                        END IF
                                    END DO
                            CLOSE(200)
    
!Constrain the randomised initial population. 
!Two DO-Loops are used to constrain the gas and steam path in HRSG. First loop initialises values at HRSG gas inlet. Second loop initialises the rest of components
                            DO I=1, NoOfPhenotypes
                                IF(StationNumber_Phenotypes(I) == 102 .AND. VariableType_Phenotypes(I) == 1) THEN !SuHe is found
                                    Phenotypes(I) = rand(Tg_exhaustGT_OD - Tg_exhaustGT_OD/NoOfHRSGComponents, Tg_exhaustGT_OD) !array index reserved for Tg_out
                                    
                                    !Search for the steam station outlet number when gas station outlet number is 102
                                    DO J=1,NoOfBricks_InputFile
                                        IF(StationNumberGas_InputFile(J,2) == 102) THEN
                                            Index2 = StationNumberSteam_InputFile(J,3)
                                        END IF
                                    END DO
                                    !Using Index2, search GA_Individual_Values until Index2 is found
                                    DO J=1, NoOfPhenotypes
                                        IF(StationNumber_Phenotypes(J) == Index2  .AND. VariableType_Phenotypes(J) == 1) THEN
                                            Phenotypes(J) = rand(Phenotypes(I), Tg_exhaustGT_OD) !array index reserved for Ts_out
                                        END IF
                                    END DO
                                    Index = I !used to transfer index to the second Do-loop
                                END IF
                            END DO

                            DO I=103, 101+NoOfHRSGComponents
                                DO J=1, NoOfPhenotypes
                                    IF(StationNumber_Phenotypes(J) == I   .AND. VariableType_Phenotypes(J) == 1) THEN
                                        Phenotypes(J) = rand(Phenotypes(Index) - Phenotypes(Index)/NoOfHRSGComponents, Phenotypes(Index)) !array index reserved for Tg_out
                                        
                                        !Search for the steam station outlet number when gas station outlet number is =I
                                        DO ii=1,NoOfBricks_InputFile
                                            IF(StationNumberGas_InputFile(ii,2) == I) THEN 
                                                Index2 = StationNumberSteam_InputFile(ii,3) !steam outlet
                                            END IF
                                        END DO
                                        !Using Index2, search GA_Individual_Values until Index2 is found
                                        DO ii=1, NoOfPhenotypes
                                            IF(StationNumber_Phenotypes(ii) == Index2   .AND. VariableType_Phenotypes(ii) == 1) THEN
                                                Phenotypes(ii) = rand(Phenotypes(J), Phenotypes(Index)) !array index reserved for Ts_out    
                                            END IF
                                        END DO
                                        Index = J
                                    END IF
                                END DO
                            END DO 
                            
                            !Constrain the Evaporator steam inlet temperature = economiser steam outlet temperature for all pressure levels
                            DO I=1,NoOfBricks_InputFile
                                IF(BrickName_InputFile(I) == "Evap") THEN
                                    DO J=1, NoOfPhenotypes
                                        IF(StationNumberSteam_InputFile(I,1) == StationNumber_Phenotypes(J) .AND. VariableType_Phenotypes(J) == 1) THEN
                                            DO ii=1, NoOfPhenotypes
                                                IF(StationNumber_Phenotypes(ii) == StationNumber_Phenotypes(J) -1 .AND. VariableType_Phenotypes(ii) == 1) THEN
                                                    Phenotypes(J) = Phenotypes(ii)
                                                END IF
                                            END DO
                                        END IF
                                    END DO
                                END IF
                            END DO
                            
!Initial constraints for pressure. Pressure levels are analysed from highest to lowest. First pressure level is a random number between 0.01 and 220. Next pressure levels are
! a random number between previous pressure and 220. Also, initialisation is done so that there is no pressure loss in the evaporator outlet.
                            Index = 0
                            DO I=101, 101+NoOfHRSGComponents
                                DO J=1, NoOfBricks_InputFile
                                    IF(StationNumberGas_InputFile(J,1) == I .AND. BrickName_InputFile(J) == "Evap") THEN
                                        DO ii=1, NoOfPhenotypes
                                            IF(StationNumber_Phenotypes(ii) == StationNumberSteam_InputFile(J,1) .AND. VariableType_Phenotypes(ii) == 2) THEN
                                                IF(Index == 0) THEN
                                                    Phenotypes(ii) = rand(0.01,220.0)
                                                    Index = ii
                                                    DO jj=1, NoOfPhenotypes
                                                        IF(StationNumber_Phenotypes(jj) == StationNumberSteam_InputFile(J,3) .AND. VariableType_Phenotypes(jj) == 2) THEN
                                                                Phenotypes(jj) = Phenotypes(ii) !initialising no pressure loss in evaporator outlet
                                                        END IF
                                                    END DO
                                                    !Search forward for SuHe brick and initialise outlet pressure = evaporator inlet pressure
                                                        Index2 = StationNumber_Phenotypes(ii)
                                                        kk = 1
                                                        DO WHILE(kk == 1)
                                                            Index2 = Index2+1
                                                            DO jj=1, NoOfPhenotypes
                                                                IF(Brick_Phenotypes(jj) == "SuHe" .AND. StationNumber_Phenotypes(jj) == Index2 .AND. VariableType_Phenotypes(jj) == 2) THEN
                                                                    Phenotypes(jj) = Phenotypes(ii)
                                                                    kk = 0
                                                                END IF
                                                            END DO
                                                        END DO
                                                    !Search backward for Econ brick and initialise outlet pressure = evaporator inlet pressure
                                                        Index2 = StationNumber_Phenotypes(ii)
                                                        kk = 1
                                                        DO WHILE(kk == 1)
                                                            Index2 = Index2-1
                                                            DO jj=1, NoOfPhenotypes
                                                                IF(Brick_Phenotypes(jj) == "Econ" .AND. StationNumber_Phenotypes(jj) == Index2 .AND. VariableType_Phenotypes(jj) == 2) THEN
                                                                    Phenotypes(jj) = Phenotypes(ii)
                                                                    kk = 0
                                                                END IF
                                                            END DO
                                                        END DO
                                                ELSE
                                                    Phenotypes(ii) = rand(0.01,Phenotypes(Index))
                                                    Index = ii
                                                    DO jj=1, NoOfPhenotypes
                                                        IF(StationNumber_Phenotypes(jj) == StationNumberSteam_InputFile(J,3) .AND. VariableType_Phenotypes(jj) == 2) THEN
                                                                Phenotypes(jj) = Phenotypes(ii) !initialising no pressure loss in evaporator
                                                        END IF
                                                    END DO
                                                    !Search forward for SuHe brick and initialise outlet pressure = evaporator inlet pressure
                                                        Index2 = StationNumber_Phenotypes(ii)
                                                        kk = 1
                                                        DO WHILE(kk == 1)
                                                            Index2 = Index2+1
                                                            DO jj=1, NoOfPhenotypes
                                                                IF(Brick_Phenotypes(jj) == "SuHe" .AND. StationNumber_Phenotypes(jj) == Index2 .AND. VariableType_Phenotypes(jj) == 2) THEN
                                                                    Phenotypes(jj) = Phenotypes(ii)
                                                                    kk = 0
                                                                END IF
                                                            END DO
                                                        END DO
                                                    !Search backward for Econ brick and initialise outlet pressure = evaporator inlet pressure
                                                        Index2 = StationNumber_Phenotypes(ii)
                                                        kk = 1
                                                        DO WHILE(kk == 1)
                                                            Index2 = Index2+1
                                                            DO jj=1, NoOfPhenotypes
                                                                IF(Brick_Phenotypes(jj) == "Econ" .AND. StationNumber_Phenotypes(jj) == Index2 .AND. VariableType_Phenotypes(jj) == 2) THEN
                                                                    Phenotypes(jj) = Phenotypes(ii)
                                                                    kk = 0
                                                                END IF
                                                            END DO
                                                        END DO
                                                END IF
                                            END IF
                                        END DO
                                    END IF
                                END DO
                            END DO



                            OPEN (UNIT=600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA_Individual_Values.txt",ACTION="WRITE")
                                WRITE(600,"(9999(G11.2,:))") "Brick", "Station_#", "Var_type","Phenotype"
                            CLOSE(600)
                            DO I=1, NoOfPhenotypes
                                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA_Individual_Values.txt",ACTION="WRITE", POSITION="APPEND")
                                    WRITE (600, "(G11.2,:)", ADVANCE='no') Brick_Phenotypes(I)
                                    WRITE (600,"(4I11)", ADVANCE='no') StationNumber_Phenotypes(I), VariableType_Phenotypes(I)
                                    WRITE (600,"(8f11.4)", ADVANCE='no') Phenotypes(I)
                                CLOSE(600)
                            END DO
 !               !Initialise the mode for db.match.CCPP_OD.exe
 !               OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\db.match.CCPP_OD\db.match.CCPP_OD\Debug\db.match.CCPP_OD_input.txt",ACTION="WRITE")
 !                           WRITE(200,*) "1 Global Mode: '0' = DP, '1' = OD(GA), '2' = OD(NR)"
 !                       CLOSE(200)
                CALL SYSTEM(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\db.match.CCPP_OD\db.match.CCPP_OD\x64\Debug\db.match.CCPP_OD.exe")
                
                OPEN (UNIT=600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA_Fitness_Function.txt",ACTION="READ")
                    READ(600,*) Fitness_PopulationArray(P)
                CLOSE(600)
    !if fitness = NaN
    IF (Fitness_PopulationArray(P) /= Fitness_PopulationArray(P)) THEN
        V = 0 !repeat initialising individual
    ELSE
        V = 1 !proceed forward
    END IF
                
    END DO !DO WHILE (V == 0)  
    
                !Write phenotypes to file
                DO I=1, NoOfPhenotypes
                    Phenotypes_PopulationArray(P,I) = Phenotypes(I)
                    Genotypes_PopulationArray(P,I) = trim(DectoBinary(INT(Phenotypes_PopulationArray(P,I)*DataResolution)))
                END DO
                
                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\1\Phenotypes.txt", ACTION="WRITE", POSITION="APPEND")
                    DO I=1, NoOfPhenotypes
                        WRITE(600,"(f11.4)", ADVANCE='no') Phenotypes_PopulationArray(P,I)
                    END DO
                CLOSE(600)
                
                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\1\Genotypes.txt", ACTION="WRITE", POSITION="APPEND")
                    DO I=1, NoOfPhenotypes
                       WRITE(600,"(A41)", ADVANCE='no') Genotypes_PopulationArray(P,I)
                    END DO
                CLOSE(600)
                
                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\1\Fitness.txt", ACTION="WRITE", POSITION="APPEND")
                       WRITE(600,"(f15.4)") Fitness_PopulationArray(P)
                CLOSE(600)     
                            !Record ff values and OD results in VARIABLES folder. 
                            OPEN (150,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_results_OD.txt",ACTION="READ")
                                READ(150,*) DummyCharacter
                                DO J=1, NoOfBricksSteamPath_results_OD
                                    READ(150,*)                         BrickName_results_OD, &
                                                                        StationNumberInlet_results_OD, &
                                                                        DummyInteger, &
                                                                        DummyInteger, &
                                                                        DummyInteger, &
                                                                        Tsin_OD, &
                                                                        Psin_OD, &
                                                                        hsin_OD, &
                                                                        Ssin_OD, &
                                                                        MFsin_OD, &
                                                                        Tsout_OD, &
                                                                        Psout_OD, &
                                                                        hsout_OD, &
                                                                        Ssout_OD, &
                                                                        MFsout_OD
                                    !check if file exists
                                    INQUIRE( file=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt", exist=dir_e )
                                    IF ( dir_e ) then
                                        !file exists
                                    ELSE
                                        OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt",ACTION="WRITE")
                                            WRITE(200,"(9999(G15.2,:))") "Tsin_OD", "Psin_OD", "hsin_OD", "Ssin_OD", "MFsin_OD", "Tsout_OD", "Psout_OD", "hsout_OD", "Ssout_OD", "MFsout_OD"
                                        CLOSE(200) 
                                    END IF
            
                                    OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt",ACTION="WRITE", POSITION="APPEND")
                                        WRITE (200,"(10f15.4)", ADVANCE='no') Tsin_OD, Psin_OD, hsin_OD, Ssin_OD, MFsin_OD, Tsout_OD, Psout_OD, hsout_OD, Ssout_OD, MFsout_OD
                                    CLOSE(200) 
        
        
                                !check if file exists. For ff record
                                    INQUIRE( file=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD_ff.txt", exist=dir_e )
                                    IF ( dir_e ) then
                                        !file exists
                                    ELSE
                                        OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD_ff.txt",ACTION="WRITE")
                                            WRITE(200,"(9999(G15.2,:))") "FF", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9"
                                        CLOSE(200) 
                                    END IF
                                    IF (BrickName_results_OD == "Econ" .OR. BrickName_results_OD == "Evap" .OR. BrickName_results_OD == "SuHe" .OR. BrickName_results_OD == "StTu") THEN
                                        OPEN (UNIT=600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA_Fitness_Function_"//trim(str(StationNumberInlet_results_OD))//".txt",ACTION="READ")
                                            READ(600,*) FitnessFunction_brick, F1, F2, F3, F4, F5, F6, F7, F8, F9
                                                OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD_ff.txt",ACTION="WRITE", POSITION="APPEND")
                                                    WRITE (200,"(9999f15.4)", ADVANCE='no') FitnessFunction_brick, F1, F2, F3, F4, F5, F6, F7, F8, F9
                                                CLOSE(200)
                                        CLOSE(600)
                                    END IF
                                END DO
        
                                !Continue to gas path    
                                READ(150,*) DummyCharacter
                                DO J=1, NoOfBricksGasPath_results_OD
                                    READ(150,*)                         BrickName_results_OD, &
                                                                        StationNumberInlet_results_OD, &
                                                                        DummyInteger, &
                                                                        Tgin_OD, &
                                                                        MFgin_OD, &
                                                                        Tgout_OD, &
                                                                        MFgout_OD
                                    !check if file exists
                                    INQUIRE( file=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt", exist=dir_e )
                                    IF ( dir_e ) then
                                        !file exists
                                    ELSE
                                        OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt",ACTION="WRITE")
                                            WRITE(200,"(9999(G15.2,:))") "Tgin_OD", "MFgin_OD", "Tgout_OD", "MFgout_OD"
                                        CLOSE(200) 
                                    END IF
                                    OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt",ACTION="WRITE", POSITION="APPEND")
                                        WRITE (200,"(5f15.4)", ADVANCE='no') Tgin_OD, MFgin_OD, Tgout_OD, MFgout_OD
                                    CLOSE(200)
                                END DO
                            CLOSE(150)
END DO !P=1, PopulationSize
        
!----------------------------------------------------- Write results to full log files
                DO P=1, PopulationSize
                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_full_phenotypes.txt", ACTION="WRITE", POSITION="APPEND")
                        DO J=1, NoOfPhenotypes
                            WRITE(600,"(f11.4)", ADVANCE='no') Phenotypes_PopulationArray(P,J)
                        END DO
                    CLOSE(600)
                END DO
                
                DO P=1, PopulationSize
                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_full_genotypes.txt", ACTION="WRITE", POSITION="APPEND")
                
                    DO J=1, NoOfPhenotypes
                       WRITE(600,"(A41)", ADVANCE='no') Genotypes_PopulationArray(P,J)
                    END DO
                    CLOSE(600)
                END DO
                
                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_full_fitness.txt", ACTION="WRITE", POSITION="APPEND")
                    DO P=1, PopulationSize
                        WRITE(600,"(f15.4)", ADVANCE='no') Fitness_PopulationArray(P)
                    END DO
                CLOSE(600)  
        

!----------------------------------------------------- END: initial population
!----------------------------------------------------- METHOD: Stochastic Selection

DO I=1, NoOfGenerations

            PRINT*, "Generation: ", I
                    !Read phenotypes, genothypes and fitness value
                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I))//"\Phenotypes.txt", ACTION="READ")
                        DO P=1, PopulationSize
                            DO J=1, NoOfPhenotypes
                                READ(600,"(f11.4)", ADVANCE='no') Phenotypes_PopulationArray(P,J)
                            END DO
                            READ(600,*) !to move the cursor to the next line
                        END DO
                    CLOSE(600)
                    
                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I))//"\Genotypes.txt", ACTION="READ")
                        DO P=1, PopulationSize
                            DO J=1, NoOfPhenotypes
                                READ(600,"(A41)", ADVANCE='no') Genotypes_PopulationArray(P,J)
                            END DO
                            READ(600,*) !to move the cursor to the next line
                        END DO
                    CLOSE(600)
                    
                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I))//"\Fitness.txt", ACTION="READ")
                        DO P=1, PopulationSize
                            READ(600,"(f15.4)") Fitness_PopulationArray(P)
                        END DO
                    CLOSE(600) 
                    
                    VariableRealDummy = ABS(MINVAL(Fitness_PopulationArray))
                    DO P=1, PopulationSize
                            !Fitness_PopulationArray(P) = ABS(Fitness_PopulationArray(P) + VariableRealDummy) !this is for maximization problem. Shift fitness array into first quadrant 
                            Fitness_PopulationArray(P) = Fitness_PopulationArray(P) ! this is for minimization problem
                    END DO
                    
                    !Sort Fitness_PopulationArray_Sorted in ascending order
                    Fitness_PopulationArray_Sorted = Fitness_PopulationArray
                    DO J=1, PopulationSize-1
                        DO P=1, PopulationSize-1
                            IF(Fitness_PopulationArray_Sorted(P) > Fitness_PopulationArray_Sorted(P+1)) THEN
                                VariableRealDummy = Fitness_PopulationArray_Sorted(P)
                                Fitness_PopulationArray_Sorted(P) = Fitness_PopulationArray_Sorted(P+1)
                                Fitness_PopulationArray_Sorted(P+1) = VariableRealDummy
                            END IF
                        END DO
                    END DO
                    !Calculate average fitness
                    VariableRealDummy = 0.0
                    DO P=1, PopulationSize
                        VariableRealDummy = VariableRealDummy + Fitness_PopulationArray(P) !to calculate sumation of all array elements
                        Average_fitness_PopulationArray = VariableRealDummy/PopulationSize
                    END DO
                    
                    !Calculate median for the fitness population
                    DummyInteger=PopulationSize/2
                    Median_fitness_PopulationArray = (Fitness_PopulationArray_Sorted(DummyInteger) + Fitness_PopulationArray_Sorted(DummyInteger+1))/2
                    
                    !calculate standard deviation
                    VariableRealDummy = 0.0
                    DO P=1, PopulationSize
                        VariableRealDummy = VariableRealDummy + (Fitness_PopulationArray_Sorted(P) - Median_fitness_PopulationArray)**2
                    END DO
                    Fitness_PopulationArray_StandardDeviation = SQRT(1/(REAL(PopulationSize)-1.0)*VariableRealDummy)
                    
                    !calculate median absolute deviation
                    VariableRealDummy = 0.0
                    DO P=1, PopulationSize
                        Fitness_PopulationArray_Sorted2(P) = ABS(Fitness_PopulationArray_Sorted(P)-Median_fitness_PopulationArray)
                    END DO
                                !Sort Fitness_PopulationArray_Sorted2 in ascending order
                                DO J=1, PopulationSize-1
                                    DO P=1, PopulationSize-1
                                        IF(Fitness_PopulationArray_Sorted2(P) > Fitness_PopulationArray_Sorted2(P+1)) THEN
                                            VariableRealDummy = Fitness_PopulationArray_Sorted2(P)
                                            Fitness_PopulationArray_Sorted2(P) = Fitness_PopulationArray_Sorted2(P+1)
                                            Fitness_PopulationArray_Sorted2(P+1) = VariableRealDummy
                                        END IF
                                    END DO
                                END DO
                    DummyInteger=PopulationSize/2
                    Fitness_PopulationArray_MedianAbsoluteDeviation = (Fitness_PopulationArray_Sorted2(DummyInteger) + Fitness_PopulationArray_Sorted2(DummyInteger+1))/2
                    !consistency constant for normal distribution. All data sets that are greater than 2 standard deviations from the mean are considered as outliners.
                    Fitness_PopulationArray_MedianAbsoluteDeviation = Fitness_PopulationArray_MedianAbsoluteDeviation*1.4826 
                    
                    DO P=1, PopulationSize
                        Fitness_PopulationArray_Sorted2(P) = ABS(Fitness_PopulationArray_Sorted(P)-Median_fitness_PopulationArray)/Fitness_PopulationArray_MedianAbsoluteDeviation
                    END DO
                                !Sort Fitness_PopulationArray_Sorted2 in ascending order
                                DO J=1, PopulationSize-1
                                    DO P=1, PopulationSize-1
                                        IF(Fitness_PopulationArray_Sorted2(P) > Fitness_PopulationArray_Sorted2(P+1)) THEN
                                            VariableRealDummy = Fitness_PopulationArray_Sorted2(P)
                                            Fitness_PopulationArray_Sorted2(P) = Fitness_PopulationArray_Sorted2(P+1)
                                            Fitness_PopulationArray_Sorted2(P+1) = VariableRealDummy
                                        END IF
                                    END DO
                                END DO
                    
                    !Calculate Q-Q plot
                    StepSize = 1/REAL(PopulationSize)
                    J = 1
                    DO P=1, PopulationSize
                        DO WHILE(ZscoresArray(J,2) <= P*StepSize .AND. J <= 818)
                            TheoreticalQuantiles(P) = ZscoresArray(J,1)
                            J = J+1
                        END DO
                    END DO
                                
                                
                    !Calculate Z-Scores. A measure of how many standard deviations is the data point from the mean.
                    !positive Z-score indicates data point is larger than mean, negative Z-score indicates it is less than mean
                    DO P=1, PopulationSize
                        SampleQuantiles(P) = (Fitness_PopulationArray_Sorted(P)-Median_fitness_PopulationArray)/Fitness_PopulationArray_StandardDeviation
                    END DO
                    !Zscore for specified X
                    Zscore_ConvergenceCriteria = (0.2-Median_fitness_PopulationArray)/Fitness_PopulationArray_StandardDeviation
                    
                    J = 1
                    DO WHILE (Zscore_ConvergenceCriteria >= ZscoresArray(J,1) .AND. J <= 818)
                        Zscore_ConvergenceCriteria_Probability = ZscoresArray(J,2)
                        J = J+1
                    END DO
                    !PRINT*, "probability is: ", Zscore_ConvergenceCriteria_Probability*100
                    
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_statistics.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600,*) "*******************************************************************************"
                                                                                                        WRITE(600,*) "*******************************************************************************"
                                                                                                        WRITE(600,*) "Generation ", I
                                                                                                        WRITE(600, *) "Mean: ", Average_fitness_PopulationArray
                                                                                                        WRITE(600, *) "Median: ", Median_fitness_PopulationArray
                                                                                                        WRITE(600, *) "Standard deviation from median: ", Fitness_PopulationArray_StandardDeviation
                                                                                                        WRITE(600, *) "Probability of fitness being below 0.2: ", Zscore_ConvergenceCriteria_Probability*100
                                                                                                        WRITE(600, *) "Median Absolute Deviation: ", Fitness_PopulationArray_MedianAbsoluteDeviation
                                                                                                    CLOSE(600)
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_statistics.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600, *) "Data"
                                                                                                        WRITE(600, *) "  -Fitness population"
                                                                                                        WRITE(600, *) "  -Median absolute deviation population"
                                                                                                        WRITE(600, *) "  -Theoretical quantiles"
                                                                                                        WRITE(600, *) "  -Sample quantiles"
                                                                                                    CLOSE(600)
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_statistics.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        DO P=1, PopulationSize    
                                                                                                            WRITE(600,"(f16.4)", ADVANCE='no') Fitness_PopulationArray_Sorted(P)
                                                                                                        END DO
                                                                                                    CLOSE(600)
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_statistics.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        DO P=1, PopulationSize    
                                                                                                            WRITE(600,"(f16.4)", ADVANCE='no') Fitness_PopulationArray_Sorted2(P)
                                                                                                        END DO
                                                                                                    CLOSE(600)
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_statistics.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        DO P=1, PopulationSize    
                                                                                                            WRITE(600,"(f16.4)", ADVANCE='no') TheoreticalQuantiles(P)
                                                                                                        END DO
                                                                                                    CLOSE(600)
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_statistics.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        DO P=1, PopulationSize    
                                                                                                            WRITE(600,"(f16.4)", ADVANCE='no') SampleQuantiles(P)
                                                                                                        END DO
                                                                                                    CLOSE(600)
                                                                                                    
                    !Add automatic winners (fi/f >1)
                    Selected = 0
                    Selection_Index = 0
                    
                    DO P=1, PopulationSize
                        !Presence = ABS(FLOOR(Fitness_PopulationArray(P)/Average_fitness_PopulationArray)) ! this is for maximization problem
                        Presence = ABS(FLOOR(Median_fitness_PopulationArray/Fitness_PopulationArray(P)))  ! this is for minimization problem
                        !Remainder(P) = (Fitness_PopulationArray(P)/Average_fitness_PopulationArray) - Presence ! this is for maximization problem
                        Remainder(P) = (Median_fitness_PopulationArray/Fitness_PopulationArray(P)) ! this is for minimization problem
                        DO J=1, Presence
                            Selected = Selected + 1
                            IF(Selected > PopulationSize) exit !for safety, when last individual is more than twice smaller than average value - presence will be more than 1 and then the Selection_Index array may go beyond its size
                            Selection_Index(Selected) = P
                            
                        END DO 
                    END DO
                                        
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600,*) "*******************************************************************************"
                                                                                                        WRITE(600,*) "*******************************************************************************"
                                                                                                        WRITE(600,*) "Generation ", I
                                                                                                        WRITE(600, *) "Population Fitness Array"
                                                                                                    CLOSE(600)
                                                                                                        DO P=1, PopulationSize
                                                                                                            OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                                WRITE(600,*) Fitness_PopulationArray(P)
                                                                                                            CLOSE(600)
                                                                                                        END DO
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600, *) "Median Fitness: ", Median_fitness_PopulationArray
                                                                                                        WRITE(600, *) "Elites: "
                                                                                                        WRITE(600,*) Selection_Index
                                                                                                        WRITE(600, *) "Selected: ", Selected
                                                                                                        WRITE(600, *) "Remainders: "
                                                                                                    CLOSE(600)
                                                                                                    DO P=1, PopulationSize
                                                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")    
                                                                                                            WRITE(600,*) Remainder(P)
                                                                                                        CLOSE(600)
                                                                                                    END DO
                    
                    !Selecting individuals used in roulette
                    !Calculate remainders in the population that were not selected as elites
                    DO P=1, PopulationSize
                        Remainder_Cumulative(P) = 0.0
                        DO J=1, P
                            Remainder_Cumulative(P) = Remainder_Cumulative(P) + Remainder(J)
                        END DO
                    END DO
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600,*) "Cumulative Remainders: "
                                                                                                    DO P=1, PopulationSize
                                                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                            WRITE(600,*) Remainder_Cumulative(P)
                                                                                                        CLOSE(600)
                                                                                                    END DO
                 
                    !Roulette wheel for remainding individuals in the population that are not elites
                    DO P=1, PopulationSize-Selected
                        !Create a random number between 0 and cumulative remainder
                        DummyReal = Remainder_Cumulative(PopulationSize)
                        RandomNumber = rand(0.0,DummyReal)
                        !increment J until correct individual is found
                        DO J=1, size(Remainder_Cumulative)
                            IF(Remainder_Cumulative(J) > RandomNumber) exit
                        END DO
                        Selection_Index(PopulationSize+1-P) = J
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600,*) "Random number: ", RandomNumber, " Index for individual selected in roulette: ", J
                                                                                                    CLOSE(600)
                    END DO
                    
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600,*) "Selection complete: "
                                                                                                        WRITE(600,*) Selection_Index
                                                                                                    CLOSE(600)
!----------------------------------------------------- END: SELECTION
!----------------------------------------------------- METHOD: Breed population to create new population (cross over)                  
                    !Shuffle selection
                    DO P=1, PopulationSize
                        J = INT(FLOOR(rand(1.0, REAL(PopulationSize))))
                        DummyInteger = Selection_Index(P)
                        Selection_Index(P) = Selection_Index(J)
                        Selection_Index(J) = DummyInteger
                    END DO
                    Phenotypes_PopulationArray_Dummy = Phenotypes_PopulationArray
                    Genotypes_PopulationArray_Dummy =Genotypes_PopulationArray
                    
                    DO P=1, PopulationSize
                        DO J=1, NoOfPhenotypes
                            Phenotypes_PopulationArray(P,J) = Phenotypes_PopulationArray_Dummy(Selection_Index(P),J)
                            Genotypes_PopulationArray(P,J) = Genotypes_PopulationArray_Dummy(Selection_Index(P),J)
                        END DO
                    END DO
                    
                    
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600,*) "--------------------- BREED"
                                                                                                        WRITE(600,*) "Shuffled selection: "
                                                                                                        WRITE(600,*) Selection_Index
                                                                                                        WRITE(600,*) "Genotypes"
                                                                                                    CLOSE(600)
                                                                                                    
                                                                                                    DO P=1, PopulationSize
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        DO J=1, NoOfPhenotypes
                                                                                                            WRITE(600,"(A41)", ADVANCE='no') Genotypes_PopulationArray(P,J)
                                                                                                        END DO
                                                                                                    CLOSE(600)
                                                                                                    END DO
                                                                                                    
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600,*) "Phenotypes"
                                                                                                    CLOSE(600)
                                                                                                    DO P=1, PopulationSize
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        DO J=1, NoOfPhenotypes
                                                                                                            WRITE(600,"(f11.4)", ADVANCE='no') Phenotypes_PopulationArray(P,J)
                                                                                                        END DO
                                                                                                    CLOSE(600)
                                                                                                    END DO
                    
                                                                            
                    !Count number of chromosomes for each genothype                                                                                 
                    DO P=1, PopulationSize
                        DO J=1, NoOfPhenotypes
                            N=0
                            DO K=1,40
                                IF(Genotypes_PopulationArray(P,J)(K:K) == "1" .AND. N == 0) THEN
                                    ChromosomeLengths(P,J) = 41-K
                                    N = 1
                                    exit
                                END IF
                            END DO      
                        END DO
                    END DO
                    
                    !Cumulative sumation of chromosomes
                    ChromosomeLengths_Cumulative = 0
                    DO P=1, PopulationSize
                        DO J=1, NoOfPhenotypes
                            DO K=1, J
                                ChromosomeLengths_Cumulative(P,J) = ChromosomeLengths_Cumulative(P,J) + ChromosomeLengths(P,K)
                            END DO
                        END DO
                    END DO
                    
    DO P=1, PopulationSize/2
        Print*, "Individuals: ", P*2-1, "&", P*2
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600,*) "*******************************************************************************"
                                                                                                        WRITE(600,*) "Individuals: ", P*2-1, "&", P*2
                                                                                                    CLOSE(600)
        condition_breed = .FALSE. !initiate condition_breed
        DO WHILE(condition_breed == .FALSE.)
                        DO J=1, NoOfPhenotypes
                            Genotypes_PopulationArray_CrossedOver(P*2-1,J) = Genotypes_PopulationArray(P*2-1,J) !for parent1
                            Genotypes_PopulationArray_CrossedOver(P*2,J) = Genotypes_PopulationArray(P*2,J) !for parent2
                        END DO
                        
                        DO J=1, NoOfPhenotypes           
                        
                                        Parent1_ChromosomeLength = ChromosomeLengths(P*2-1,J)
                                        Parent2_ChromosomeLength = ChromosomeLengths(P*2,J)
                                        
                                        !DummyInteger = the longer chomosome length
                                        IF (Parent1_ChromosomeLength > Parent2_ChromosomeLength) THEN
                                            DummyInteger = Parent1_ChromosomeLength
                                        ELSE
                                            DummyInteger = Parent2_ChromosomeLength
                                        END IF


                                        !Chromosome lengths for partnet1 and parent2 may be of different lengths. 
                                        !To ensure that the random number for the split is not larger than the shorter chromosome length
                                        CrossOver_Split2 = 99999 !arbitrary large number to triger do while loop below 
                                        DO WHILE(CrossOver_Split2 > Parent1_ChromosomeLength .AND. CrossOver_Split2 > Parent2_ChromosomeLength)
                                            CrossOver_Split1 = rand(1.0, REAL(DummyInteger))
                                            CrossOver_Split2 = rand(1.0, REAL(DummyInteger))
                        
                                            !sort so that Split1 < Split2
                                            IF (CrossOver_Split1 > CrossOver_Split2) THEN
                                                K = CrossOver_Split1
                                                CrossOver_Split1 = CrossOver_Split2
                                                CrossOver_Split2 = K
                                            END IF
                                        END DO
                                                                                                                   OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                                        WRITE(600,*) "Phenotype ", J, " Split points for annulus: ", CrossOver_Split1, " | ", CrossOver_Split2
                                                                                                                    CLOSE(600)
                                                                                                                    
                                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                                        WRITE(600,*) " ----- Old1 Old2 ", BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2-1,J))/DataResolution, BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2,J))/DataResolution
                                                                                                                    CLOSE(600)        
                                        DO K=40-DummyInteger+CrossOver_Split2, 40-DummyInteger+CrossOver_Split1, -1
                                            C_Parent1 = Genotypes_PopulationArray_CrossedOver(P*2-1,J)(K:K)
                                            C_Parent2 = Genotypes_PopulationArray_CrossedOver(P*2,J)(K:K)
                                            Genotypes_PopulationArray_CrossedOver(P*2-1,J)(K:K) = C_Parent2
                                            Genotypes_PopulationArray_CrossedOver(P*2,J)(K:K) = C_Parent1
                                        END DO
                                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                                        WRITE(600,*) " ----- New1 New2 ", BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2-1,J))/DataResolution, BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2,J))/DataResolution
                                                                                                                    CLOSE(600)  
                        
                        END DO !number of phenotypes
                                                                                                    !DO K=1, PopulationSize
                                                                                                    !OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                    !    DO J=1, NoOfPhenotypes
                                                                                                    !        WRITE(600,"(A41)", ADVANCE='no') Genotypes_PopulationArray(K,J)
                                                                                                    !    END DO
                                                                                                    !CLOSE(600)
                                                                                                    !END DO
                        
                                                                                                    !DO K=1, PopulationSize
                                                                                                    !OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                    !    DO J=1, NoOfPhenotypes
                                                                                                    !        Phenotypes_PopulationArray(K,J) = BinarytoDec(Genotypes_PopulationArray(K,J))/DataResolution
                                                                                                    !        WRITE(600,"(f11.4)", ADVANCE='no') Phenotypes_PopulationArray(K,J)
                                                                                                    !    END DO
                                                                                                    !CLOSE(600)
                                                                                                    !END DO
                        
!----------------------------------------------------- END: CROSS OVER 
!----------------------------------------------------- METHOD: Mutate genotypes to diversify the population
                        
                        DO Index=1, NoOfPhenotypes
                            !mutating for Child1
                            DO J=40-ChromosomeLengths(P*2-1,Index), 40
                                RandomNumber = rand(0.0,1.0)

                                IF(RandomNumber < ProbOfMutation .AND. Genotypes_PopulationArray_CrossedOver(P*2-1,Index)(J:J) == "0") THEN
                                                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                            WRITE(600,*) "Mutating individual ",  P*2-1, "genotype ", Index
                                                                                                            WRITE(600,*) "Old genotype ", trim(Genotypes_PopulationArray_CrossedOver(P*2-1,Index)), " | ", BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2-1,Index))/DataResolution
                                                                                                        CLOSE(600)
                                        Genotypes_PopulationArray_CrossedOver(P*2-1,Index)(J:J) = "1"
                                                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                            WRITE(600,*) "Gene number ", J, " replaced by new gene ", Genotypes_PopulationArray_CrossedOver(P*2-1,Index)(J:J)
                                                                                                            WRITE(600,*) "New genotype ", trim(Genotypes_PopulationArray_CrossedOver(P*2-1,Index)), " | ", BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2-1,Index))/DataResolution
                                                                                                        CLOSE(600)
                                ELSE IF(RandomNumber < ProbOfMutation .AND. Genotypes_PopulationArray_CrossedOver(P*2-1,Index)(J:J) == "1") THEN
                                                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                            WRITE(600,*) "Mutating individual ",  P*2-1, "genotype ", Index
                                                                                                            WRITE(600,*) "Old genotype ", trim(Genotypes_PopulationArray_CrossedOver(P*2-1,Index)), " | ", BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2-1,Index))/DataResolution
                                                                                                        CLOSE(600)
                                        Genotypes_PopulationArray_CrossedOver(P*2-1,Index)(J:J) = "0"
                                                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                            WRITE(600,*) "Gene number ", J, " replaced by new gene ", Genotypes_PopulationArray_CrossedOver(P*2-1,Index)(J:J)
                                                                                                            WRITE(600,*) "New genotype ", trim(Genotypes_PopulationArray_CrossedOver(P*2-1,Index)), " | ", BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2-1,Index))/DataResolution
                                                                                                        CLOSE(600)
                                END IF    
                            END DO
                            !mutating for Child2    
                            DO J=40-ChromosomeLengths(P*2,Index), 40
                                RandomNumber = rand(0.0,1.0)
                                IF(RandomNumber < ProbOfMutation .AND. Genotypes_PopulationArray_CrossedOver(P*2,Index)(J:J) == "0") THEN
                                                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                            WRITE(600,*) "Mutating individual ",  P*2, "genotype ", Index
                                                                                                            WRITE(600,*) "Old genotype ", trim(Genotypes_PopulationArray_CrossedOver(P*2,Index)), " | ", BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2,Index))/DataResolution
                                                                                                        CLOSE(600)
                                        Genotypes_PopulationArray_CrossedOver(P*2,Index)(J:J) = "1"
                                                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                            WRITE(600,*) "Gene number ", J, " replaced by new gene ", Genotypes_PopulationArray_CrossedOver(P*2,Index)(J:J)
                                                                                                            WRITE(600,*) "New genotype ", trim(Genotypes_PopulationArray_CrossedOver(P*2,Index)), " | ", BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2,Index))/DataResolution
                                                                                                        CLOSE(600)
                                ELSE IF(RandomNumber < ProbOfMutation .AND. Genotypes_PopulationArray_CrossedOver(P*2,Index)(J:J) == "1") THEN
                                                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                            WRITE(600,*) "Mutating individual ",  P*2, "genotype ", Index
                                                                                                            WRITE(600,*) "Old genotype ", trim(Genotypes_PopulationArray_CrossedOver(P*2,Index)), " | ", BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2,Index))/DataResolution
                                                                                                        CLOSE(600)
                                        Genotypes_PopulationArray_CrossedOver(P*2,Index)(J:J) = "0"
                                                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                            WRITE(600,*) "Gene number ", J, " replaced by new gene ", Genotypes_PopulationArray_CrossedOver(P*2,Index)(J:J)
                                                                                                            WRITE(600,*) "New genotype ", trim(Genotypes_PopulationArray_CrossedOver(P*2,Index)), " | ", BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2,Index))/DataResolution
                                                                                                        CLOSE(600)
                                END IF    
                            END DO
                        END DO
                            

                        
!----------------------------------------------------- END: MUTATION

                        !convert genotypes to phenotypes
                        DO J=1, NoOfPhenotypes
                            Phenotypes_PopulationArray_CrossedOver(P*2-1,J) = BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2-1,J))/DataResolution !for child1
                            Phenotypes_PopulationArray_CrossedOver(P*2,J) = BinarytoDec(Genotypes_PopulationArray_CrossedOver(P*2,J))/DataResolution !for child2
                        END DO

!----------------------------------------------------- CHECK CONSTRAINTS for newely generated phenotypes
                    DO S=1, 0, -1 !decreasing loop starts with repeat for S=1 for child1 then S=0 for child 2   
                        
                        !If Phenotype = 0, then db.match.CCPP_OD.exe will not converge. Therefore a hard constrain is given here.
                        DO J=1, NoOfPhenotypes
                            IF (Phenotypes_PopulationArray_CrossedOver(P*2-S,J) == 0) THEN
                                Phenotypes_PopulationArray_CrossedOver(P*2-S,J) = 1
                                Genotypes_PopulationArray_CrossedOver(P*2-S,J) = trim(DectoBinary(INT(Phenotypes_PopulationArray_CrossedOver(P*2-S,J)*DataResolution)))
                            END IF
                        END DO
                        
                    
!----------------------------------------------------- END: CHECK CONSTRAINTS for newely generated phenotypes        




                        OPEN (UNIT=600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA_Individual_Values.txt",ACTION="WRITE")
                            WRITE(600,"(9999(G11.2,:))") "Brick", "Station_#", "Var_type","Phenotype"
                        CLOSE(600)
                                DO J=1, NoOfPhenotypes
                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA_Individual_Values.txt",ACTION="WRITE", POSITION="APPEND")
                                        WRITE (600, "(G11.2,:)", ADVANCE='no') Brick_Phenotypes(J)
                                        WRITE (600,"(4I11)", ADVANCE='no') StationNumber_Phenotypes(J), VariableType_Phenotypes(J)
                                        WRITE (600,"(8f11.4)", ADVANCE='no') Phenotypes_PopulationArray_CrossedOver(P*2-S,J)
                                    CLOSE(600)
                                END DO
                        
 

                        
!-----------------------------------------------
!---------- interupt results of GA, and run NR.
!-----------------------------------------------
        OPEN (UNIT=600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\NR_Individual_Values.txt",ACTION="WRITE")
            WRITE(600,"(9999(G11.2,:))") "Brick", "Station_#", "Var_type","Phenotype"
        CLOSE(600)
                DO J=1, NoOfPhenotypes
                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\NR_Individual_Values.txt",ACTION="WRITE", POSITION="APPEND")
                        WRITE (600, "(G11.2,:)", ADVANCE='no') Brick_Phenotypes(J)
                        WRITE (600,"(4I11)", ADVANCE='no') StationNumber_Phenotypes(J), VariableType_Phenotypes(J)
                        WRITE (600,"(8f11.4)", ADVANCE='no') Phenotypes_PopulationArray_CrossedOver(P*2-S,J)
                    CLOSE(600)
                END DO
        !CALL SYSTEM(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\NR\NR\Debug\NR.exe")
!-----------------------------------------------
!---------- end NR.
!-----------------------------------------------

                        !Initialise the mode for db.match.CCPP_OD.exe (in order to continue running GA)
                        OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\db.match.CCPP_OD\db.match.CCPP_OD\x64\Debug\db.match.CCPP_OD_input.txt",ACTION="WRITE")
                            WRITE(200,*) "1 Global Mode: '0' = DP, '1' = OD(GA), '2' = OD(NR)"
                        CLOSE(200)  
                        
                        CALL SYSTEM(trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\db.match.CCPP_OD\db.match.CCPP_OD\x64\Debug\db.match.CCPP_OD.exe")
                        
                        OPEN (UNIT=600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA_Fitness_Function.txt",ACTION="READ")
                            READ(600,*) Fitness_PopulationArray(P*2-S)
                        CLOSE(600)
                        
                        condition_breed = .TRUE. !initiate .TRUE. value now. The following code can be used to introduce constraints. if .FALSE. the GA will generate new breed for Child1 and Child2.
                        !if fitness = NaN
                        IF(Fitness_PopulationArray(P*2-S) /= Fitness_PopulationArray(P*2-S)) THEN 
                            print*, "Fitness is NaN"
                            Fitness_PopulationArray(P*2-S) = 9999
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600,*) "ERROR - fitness is NaN. -----------------------"
                                                                                                    CLOSE(600) 
                            !condition_breed = .FALSE.
                        END IF
                        
                        !if child1 or child2 fitness > 1000. Sometimes the fitness is very large, then a number doesnt fit in the pre-located space in text file - showing ********* instead of the number itself
                        IF(Fitness_PopulationArray(P*2-S) > 9999) THEN
                            print*, "Fitness is above 9,999"
                            Fitness_PopulationArray(P*2-S) = 9999
                                                                                                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                                                                                                        WRITE(600,*) "ERROR - fitness is > 9,999 ---------------------"
                                                                                                    CLOSE(600)
                            !condition_breed = .FALSE.
                        END IF
                        
                        
                        
                        
                        !if child1 or child2 fitness > 1000. Sometimes the fitness is very large, then a number doesnt fit in the pre-located space in text file - showing ********* instead of the number itself
                        !IF(Fitness_PopulationArray(P*2-1) > 9999 .OR. Fitness_PopulationArray(P*2) > 9999) THEN
                        !    print*, "Fitness is above 9,999"
                        !                                                                            OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                        !                                                                                WRITE(600,*) "ERROR - fitness is > 9,999 ---------------------"
                         !                                                                           CLOSE(600)
                            !condition_breed = .FALSE.
                        !END IF
                        
                        
                        !if child1 or child2 fitness = NaN
                        !IF(Fitness_PopulationArray(P*2-1) /= Fitness_PopulationArray(P*2-1) .OR. Fitness_PopulationArray(P*2) /= Fitness_PopulationArray(P*2)) THEN 
                            !PRINT*, "Fitness is NaN. Generation ", I, " Individuals ", P*2-1, " & ", P*2
                            !PRINT*, "Child1 ", Fitness_PopulationArray(P*2-1)
                            !PRINT*, "Child2 ", Fitness_PopulationArray(P*2)
                            !PRINT*, "Child1 phenotype "
                            !DO J=1, NoOfPhenotypes
                            !    PRINT*, Phenotypes_PopulationArray_CrossedOver(P*2-1,J)
                            !END DO
                            !PRINT*, "Child2 phenotype "
                            !DO J=1, NoOfPhenotypes
                            !    PRINT*, Phenotypes_PopulationArray_CrossedOver(P*2,J)
                            !END DO
                            !PRINT*, "Fitness is NaN. Trying again"
                            !condition_breed = .FALSE.
                            !                                                                        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
                            !                                                                            WRITE(600,*) "ERROR - fitness is NaN. -----------------------"
                            !                                                                        CLOSE(600) 
                            !pause
                        !ELSE
                        !    condition_breed = .TRUE.
                        !END IF
                        
                        
                        
                        !If child1 and child2 have passed all constraints above, condition_breed will be .TRUE.
                        !only then record values in logs
                        IF ( condition_breed == .TRUE. ) THEN
                            OPEN (150,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_results_OD.txt",ACTION="READ")
                                READ(150,*) DummyCharacter
                                DO J=1, NoOfBricksSteamPath_results_OD
                                    READ(150,*)                         BrickName_results_OD, &
                                                                        StationNumberInlet_results_OD, &
                                                                        DummyInteger, &
                                                                        DummyInteger, &
                                                                        DummyInteger, &
                                                                        Tsin_OD, &
                                                                        Psin_OD, &
                                                                        hsin_OD, &
                                                                        Ssin_OD, &
                                                                        MFsin_OD, &
                                                                        Tsout_OD, &
                                                                        Psout_OD, &
                                                                        hsout_OD, &
                                                                        Ssout_OD, &
                                                                        MFsout_OD
                                    !check if file exists
                                    INQUIRE( file=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt", exist=dir_e )
                                    IF ( dir_e ) then
                                        !file exists
                                    ELSE
                                        OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt",ACTION="WRITE")
                                            WRITE(200,"(9999(G15.2,:))") "Tsin_OD", "Psin_OD", "hsin_OD", "Ssin_OD", "MFsin_OD", "Tsout_OD", "Psout_OD", "hsout_OD", "Ssout_OD", "MFsout_OD"
                                        CLOSE(200) 
                                    END IF
            
                                    OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt",ACTION="WRITE", POSITION="APPEND")
                                        WRITE (200,"(10f15.4)", ADVANCE='no') Tsin_OD, Psin_OD, hsin_OD, Ssin_OD, MFsin_OD, Tsout_OD, Psout_OD, hsout_OD, Ssout_OD, MFsout_OD
                                    CLOSE(200) 
        
        
                                !check if file exists. For ff record
                                    INQUIRE( file=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD_ff.txt", exist=dir_e )
                                    IF ( dir_e ) then
                                        !file exists
                                    ELSE
                                        OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD_ff.txt",ACTION="WRITE")
                                            WRITE(200,"(9999(G15.2,:))") "FF", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9"
                                        CLOSE(200) 
                                    END IF
                                    IF (BrickName_results_OD == "Econ" .OR. BrickName_results_OD == "Evap" .OR. BrickName_results_OD == "SuHe" .OR. BrickName_results_OD == "StTu") THEN
                                        OPEN (UNIT=600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA_Fitness_Function_"//trim(str(StationNumberInlet_results_OD))//".txt",ACTION="READ")
                                            READ(600,*) FitnessFunction_brick, F1, F2, F3, F4, F5, F6, F7, F8, F9
                                                OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD_ff.txt",ACTION="WRITE", POSITION="APPEND")
                                                    WRITE (200,"(9999f15.4)", ADVANCE='no') FitnessFunction_brick, F1, F2, F3, F4, F5, F6, F7, F8, F9
                                                CLOSE(200)
                                        CLOSE(600)
                                    END IF
                                END DO
        
                                !Continue to gas path    
                                READ(150,*) DummyCharacter
                                DO J=1, NoOfBricksGasPath_results_OD
                                    READ(150,*)                         BrickName_results_OD, &
                                                                        StationNumberInlet_results_OD, &
                                                                        DummyInteger, &
                                                                        Tgin_OD, &
                                                                        MFgin_OD, &
                                                                        Tgout_OD, &
                                                                        MFgout_OD
                                    !check if file exists
                                    INQUIRE( file=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt", exist=dir_e )
                                    IF ( dir_e ) then
                                        !file exists
                                    ELSE
                                        OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt",ACTION="WRITE")
                                            WRITE(200,"(9999(G15.2,:))") "Tgin_OD", "MFgin_OD", "Tgout_OD", "MFgout_OD"
                                        CLOSE(200) 
                                    END IF
                                    OPEN (200,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Variables\"//trim(str(StationNumberInlet_results_OD))//"_results_OD.txt",ACTION="WRITE", POSITION="APPEND")
                                        WRITE (200,"(5f15.4)", ADVANCE='no') Tgin_OD, MFgin_OD, Tgout_OD, MFgout_OD
                                    CLOSE(200)
                                END DO
                            CLOSE(150)
                        END IF
                       
                        
!if this pause is engaged, one can monitor every individual in the population using MATLAB program to see changes in values of station variables across population
!pause 
                   
              END DO !end S loop. S=0 for child1 and S=1 for child 2        

        END DO !end a condition_breed loop. Break out of loop if condition_breed = .TRUE.
                       
    END DO !end P loop for PopulationSize/2                


                    
!----------------------------------------------------- Update the LOG

    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
        WRITE(600,*) "Genotypes crossed over"
    CLOSE(600)
    DO P=1, PopulationSize
        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
            DO J=1, NoOfPhenotypes
                WRITE(600,"(A41)", ADVANCE='no') Genotypes_PopulationArray_CrossedOver(P,J)
            END DO
        CLOSE(600)
    END DO
    
    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
        WRITE(600,*) "Phenotypes crossed over"
    CLOSE(600)
    DO P=1, PopulationSize
        OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log.txt", ACTION="WRITE", POSITION="APPEND")
            DO J=1, NoOfPhenotypes
                WRITE(600,"(f11.4)", ADVANCE='no') Phenotypes_PopulationArray_CrossedOver(P,J)
            END DO
        CLOSE(600)
    END DO

             
               
!----------------------------------------------------- Write results into the next generation
                DO P=1, PopulationSize

                    
                    
                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I+1))//"\Phenotypes.txt", ACTION="WRITE", POSITION="APPEND")
                        DO J=1, NoOfPhenotypes
                            WRITE(600,"(f11.4)", ADVANCE='no') Phenotypes_PopulationArray_CrossedOver(P,J)
                        END DO
                    CLOSE(600)
                
                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I+1))//"\Genotypes.txt", ACTION="WRITE", POSITION="APPEND")
                        DO J=1, NoOfPhenotypes
                           WRITE(600,"(A41)", ADVANCE='no') Genotypes_PopulationArray_CrossedOver(P,J)
                        END DO
                    CLOSE(600)
                
                    OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\GA\Generations\"//trim(str(I+1))//"\Fitness.txt", ACTION="WRITE", POSITION="APPEND")
                           WRITE(600,"(f15.4)") Fitness_PopulationArray(P)
                    CLOSE(600)          
                    
                    
                END DO
                
                
!----------------------------------------------------- Write results to full log files
                DO P=1, PopulationSize
                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_full_phenotypes.txt", ACTION="WRITE", POSITION="APPEND")
                    DO J=1, NoOfPhenotypes
                        WRITE(600,"(f11.4)", ADVANCE='no') Phenotypes_PopulationArray_CrossedOver(P,J)
                    END DO
                CLOSE(600)
                END DO
                
                DO P=1, PopulationSize
                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_full_genotypes.txt", ACTION="WRITE", POSITION="APPEND")
                
                    DO J=1, NoOfPhenotypes
                       WRITE(600,"(A41)", ADVANCE='no') Genotypes_PopulationArray_CrossedOver(P,J)
                    END DO
                CLOSE(600)
                END DO
                
                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_full_fitness.txt", ACTION="WRITE", POSITION="APPEND")
                DO P=1, PopulationSize
                    WRITE(600,"(f15.4)", ADVANCE='no') Fitness_PopulationArray(P)
                END DO
                CLOSE(600)  
                                                                     
!-------------------------------- write MINIMUM fitness results to a file
                !Search for minimum fitness
                VariableRealDummy = 10000
                do P=1, PopulationSize
                    if(Fitness_PopulationArray(P) < VariableRealDummy) then
                        VariableRealDummy = Fitness_PopulationArray(P)
                        Index = P
                    end if
                end do
                !----------------------------------------------------- Write minimum results to full log files
                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_min_phenotypes.txt", ACTION="WRITE", POSITION="APPEND")
                    DO J=1, NoOfPhenotypes
                        WRITE(600,"(f11.4)", ADVANCE='no') Phenotypes_PopulationArray_CrossedOver(Index,J)
                    END DO
                CLOSE(600)
                
                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_min_genotypes.txt", ACTION="WRITE", POSITION="APPEND")
                    DO J=1, NoOfPhenotypes
                       WRITE(600,"(A41)", ADVANCE='no') Genotypes_PopulationArray_CrossedOver(Index,J)
                    END DO
                CLOSE(600)
                
                OPEN (600,FILE=trim(db_match_CCPP_file_path)//"\db.match.CCPP_OD\Log_min_fitness.txt", ACTION="WRITE", POSITION="APPEND")
                    WRITE(600,"(f15.4)", ADVANCE='no') Fitness_PopulationArray(Index)
                CLOSE(600)  

END DO !end I loop for NoOfGenerations
        
        
        

    
!------------------------------------------------------------- Write results
    CALL WriteVar
    
    
    
    
    
pause
    
end program MainProgram


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! PROGRAM END !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------
character(len=20) function str(k) !Function to convert an integer to string
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
    
end function str
!---------------------------------------------------------------------------------------
REAL function rand(min,max)

    real, intent(in) :: min
    real, intent(in) :: max
    real             :: randnum
    CALL random_number(randnum)
    rand = min + randnum*(max - min)
    !Print*, "random number is: ", rand
end function rand
!---------------------------------------------------------------------------------------
CHARACTER(len=40) function DectoBinary(input)
IMPLICIT NONE
!

INTEGER, INTENT(IN) :: input  !!ARGUMENT COMING IN TO FUNCTION
INTEGER(8), DIMENSION(30) :: base2_array !there is an upper bond on the amount of numbers for the integer, i.e. 2^31 and higher will be too high number for the integer to hold
INTEGER :: I, J
INTEGER, ALLOCATABLE, DIMENSION(:) :: integerArray
INTEGER, DIMENSION(40) :: output_integerArray
INTEGER :: quotient
CHARACTER(1), DIMENSION(40) :: output_stringArray
CHARACTER(len=40) :: output_string !assuming the lenght
!CHARACTER(len=20) :: str !assuming the lenght
INTEGER :: output_integer

!I am getting error when input = 2048. problem fixed

base2_array = 0
I = 0
J = 0
quotient = 0
output_string = ""
output_integer = 0
output_integerArray = 0
output_stringArray = ""

    !input = 8192
    
    DO I=1, 30
        base2_array(30-(I-1)) = 2**(I-1)
    END DO
    
    DO I=1, 30
        IF((input-base2_array(I)) >= 0) THEN
            J = I
            exit
        END IF
    END DO
    ALLOCATE (integerArray(30-J+1))
    integerArray = 0

    
    quotient = input
    J = size(integerArray)

    
    DO I=1, 30
        IF((quotient-base2_array(I)) >= 0) THEN
            integerArray(J-(30-I)) = 1
            quotient = quotient - base2_array(I)
        END IF
    END DO
    

    
    DO I=41-size(integerArray), 40
        output_integerArray(I) = integerArray(I+size(integerArray)-40)
    END DO
    
        !convert integer array to string array
    DO I=1, size(output_integerArray)
        write(output_stringArray(I), '(I1)')  output_integerArray(I)
    END DO
    
    
    !convert string array to character
    DO I=1, size(output_stringArray)
        output_string(I:I) = output_stringArray(I)
    END DO
    
    
    !write(6,*) trim(string), "bla ", size(output_integerArray)
    
    !convert integer to string
    !write (str, *) size(output_integerArray)
    !str = adjustl(str)
   
    !convert character into integer
    
    !read(string,"(I"//trim(str)//")") output_integer
    !read(string,"(I10)") output_integer
    !DO I=1, size(output_stringArray)
     !   read(output_string(1:I), '(I)') output_integer
    !END DO
    
    DectoBinary = output_string
END function DectoBinary
!---------------------------------------------------------------------------------------
REAL function BinarytoDec(input)
IMPLICIT NONE
!

CHARACTER(len=40), INTENT(IN) :: input  !!ARGUMENT COMING IN TO FUNCTION
INTEGER(8), DIMENSION(30) :: base2_array !there is an upper bond on the amount of numbers for the integer, i.e. 2^31 and higher will be too high number for the integer to hold
INTEGER :: I, J, K
INTEGER, ALLOCATABLE, DIMENSION(:) :: input_integerArray
INTEGER :: input_length
INTEGER :: output_integer 
INTEGER :: index 
base2_array = 0
I = 0
J = 0
K = 0
output_integer = 0
input_length = 0
output_integer = 0
index = 1
!getting error when input = 0000000000000000000000000000000000000000. Problem fixed


    DO I=1, 30
        base2_array(30-(I-1)) = 2**(I-1)
    END DO

    !count length of binary string
    K=0
    DO I=1, 40
        IF(input(I:I) == "1" .AND. K == 0) THEN
            input_length = 41-I
            index = I
            K = 1
            exit
        END IF
    END DO
    IF ( index == 1) THEN
        input_length = 40
    END IF
    
    ALLOCATE (input_integerArray(input_length))
    input_integerArray = 0
    !convert string into integer array
    DO I=index,40 
        read(input(I:I), '(I)') input_integerArray(I-index+1)
    END DO


    DO I=1, input_length
        IF(input_integerArray(I) == 1) THEN
            output_integer = output_integer + base2_array(30-input_length+I)
        END IF
    END DO
    
    BinarytoDec = output_integer
    
END function BinarytoDec





! split a string into 2 either side of a delimiter token
  !SUBROUTINE split_string(instring, string1, string2, delim)
  !!  CHARACTER(30) :: instring,delim
   ! CHARACTER(30),INTENT(OUT):: string1,string2
   ! INTEGER :: index

!    instring = TRIM(instring)

!    index = SCAN(instring,delim)
!    string1 = instring(1:index-1)
!    string2 = instring(index+1:)

 ! END SUBROUTINE split_string